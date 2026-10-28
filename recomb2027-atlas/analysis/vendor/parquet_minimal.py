"""Strict offline reader for the flat, SNAPPY, v1 pages in supplied Parquet crops.
Not a general Parquet implementation. Raises on unsupported schemas/encodings.
Prefer pandas+pyarrow when installed. Apache Thrift handles footer/page metadata;
the system libsnappy handles decompression. Used only for dependency recovery.
"""
from pathlib import Path
import struct, ctypes, ctypes.util
import numpy as np
import pandas as pd
from thrift.protocol import TCompactProtocol
from thrift.transport import TTransport
from thrift.Thrift import TType

def _read(p,t):
    if t==TType.STRUCT:
        d={};p.readStructBegin()
        while True:
            _,kind,n=p.readFieldBegin()
            if kind==TType.STOP:break
            d[n]=_read(p,kind);p.readFieldEnd()
        p.readStructEnd();return d
    if t in (TType.LIST,TType.SET):
        kind,n=p.readListBegin();x=[_read(p,kind) for _ in range(n)];p.readListEnd();return x
    if t==TType.MAP:
        k,v,n=p.readMapBegin();x={_read(p,k):_read(p,v) for _ in range(n)};p.readMapEnd();return x
    return {TType.BOOL:p.readBool,TType.BYTE:p.readByte,TType.I16:p.readI16,TType.I32:p.readI32,TType.I64:p.readI64,TType.DOUBLE:p.readDouble,TType.STRING:p.readBinary}[t]()

def _thrift(buf):
    t=TTransport.TMemoryBuffer(buf);d=_read(TCompactProtocol.TCompactProtocol(t),TType.STRUCT)
    return d,t._buffer.tell()

_LIB=ctypes.CDLL(ctypes.util.find_library('snappy'))
_LIB.snappy_uncompress.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_void_p,ctypes.POINTER(ctypes.c_size_t)]
_LIB.snappy_uncompress.restype=ctypes.c_int

def _unpack(data,codec,n):
    if codec==0:
        if len(data)!=n:raise ValueError('bad uncompressed size')
        return data
    if codec!=1:raise NotImplementedError(f'codec {codec}')
    out=ctypes.create_string_buffer(n);size=ctypes.c_size_t(n)
    status=_LIB.snappy_uncompress(data,len(data),out,ctypes.byref(size))
    if status or size.value!=n:raise ValueError('SNAPPY decompression failed')
    return out.raw[:n]

def _varint(data,p):
    x=0;s=0
    while True:
        b=data[p];p+=1;x|=(b&127)<<s
        if b<128:return x,p
        s+=7
        if s>63:raise ValueError('varint overflow')

def _rle(data,width,n):
    ans=[];p=0
    if not 0<=width<=64:raise ValueError('bad bit width')
    while len(ans)<n:
        header,p=_varint(data,p)
        if header==0:raise ValueError('zero run')
        if header&1:
            count=(header>>1)*8;size=count*width//8
            raw=int.from_bytes(data[p:p+size],'little');p+=size;mask=(1<<width)-1
            ans.extend((raw>>(i*width))&mask for i in range(count))
        else:
            count=header>>1;size=(width+7)//8;value=int.from_bytes(data[p:p+size],'little');p+=size
            ans.extend([value]*count)
    return ans[:n]

def _plain(data,typ,n):
    if typ==0:return np.unpackbits(np.frombuffer(data,dtype=np.uint8),bitorder='little')[:n].astype(bool).tolist()
    dtypes={1:'<i4',2:'<i8',4:'<f4',5:'<f8'}
    if typ in dtypes:
        ar=np.frombuffer(data,dtype=dtypes[typ],count=n)
        return ar.tolist()
    if typ==6:
        values=[];p=0
        for _ in range(n):
            length=struct.unpack_from('<I',data,p)[0];p+=4
            values.append(data[p:p+length].decode('utf-8'));p+=length
        if p>len(data):raise ValueError('truncated binary')
        return values
    raise NotImplementedError(f'physical type {typ}')

def read_flat_parquet(path):
    b=Path(path).read_bytes()
    if b[:4]!=b'PAR1' or b[-4:]!=b'PAR1':raise ValueError('not parquet')
    size=struct.unpack('<I',b[-8:-4])[0];meta,_=_thrift(b[-size-8:-8]);schema=meta[2]
    if len(schema)!=schema[0][5]+1:raise NotImplementedError('nested schema')
    specs={s[4].decode():s for s in schema[1:]};output={name:[] for name in specs}
    for group in meta[4]:
        for col in group[1]:
            m=col[3];name=m[3][0].decode();spec=specs[name]
            if len(m[3])!=1 or spec[3] not in (0,1):raise NotImplementedError('nested/repeated')
            count=m[5];p=m.get(11,m[9]);dictionary=None;values=[]
            while len(values)<count:
                head,k=_thrift(b[p:]);p+=k
                raw=_unpack(b[p:p+head[3]],m[4],head[2]);p+=head[3]
                if head[1]==2:
                    d=head[7]
                    if d[2]!=0:raise NotImplementedError('dictionary encoding')
                    dictionary=_plain(raw,m[1],d[1]);continue
                if head[1]!=0:raise NotImplementedError('not V1 data page')
                h=head[5];n=h[1];offset=0
                if spec[3]==1:
                    length=struct.unpack_from('<I',raw,offset)[0];offset+=4
                    defined=_rle(raw[offset:offset+length],1,n);offset+=length
                else:defined=[1]*n
                present=sum(defined);payload=raw[offset:]
                if h[2] in (2,8):
                    if dictionary is None:raise ValueError('missing dictionary')
                    indices=_rle(payload[1:],payload[0],present)
                    chunk=[dictionary[j] for j in indices]
                elif h[2]==0:chunk=_plain(payload,m[1],present)
                else:raise NotImplementedError(f'encoding {h[2]}')
                it=iter(chunk);values.extend(next(it) if flag else None for flag in defined)
            if len(values)!=count:raise ValueError('count mismatch')
            output[name].extend(values)
    if any(len(v)!=meta[3] for v in output.values()):raise ValueError('row count mismatch')
    return pd.DataFrame(output)
