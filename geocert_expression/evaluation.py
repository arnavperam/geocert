"""Common evaluation on supplied, aligned ownership truth and interval outputs.

No method is run by this module. Cell labels must already be mapped to the SAME
biological-cell identities; comparing arbitrary segmentation IDs is invalid.
Truth from a second segmenter is an annotation comparator, not anatomical truth.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from .calibration import wilson

def misassignment_auc(fragilities, errors):
    if len(fragilities)!=len(errors):raise ValueError('score/label length mismatch')
    pos=[s for s,e in zip(fragilities,errors) if e]
    neg=[s for s,e in zip(fragilities,errors) if not e]
    if not pos or not neg:return None
    return float(sum((a>b)+F(1,2)*(a==b) for a in pos for b in neg)/(len(pos)*len(neg)))

def evaluate_aligned(records, intervals):
    """RNA record keys: region,id,gene,truth_owner,predicted_owner,fragility.

    Interval keys: region,cell,gene,lower,upper. All tested count keys must be
    supplied, including zero-count truth entries represented by an upper bound.
    Coverage is reported at independent region unit, not per-transcript binomial.
    """
    truth=Counter();seen=set();errors=[];scores=[];regions=set()
    for r in records:
        key=(r['region'],r['id'])
        if key in seen:raise ValueError('duplicate transcript in evaluation region')
        seen.add(key);regions.add(r['region'])
        truth[(r['region'],r['truth_owner'],r['gene'])]+=1
        errors.append(r['truth_owner']!=r['predicted_owner']);scores.append(r['fragility'])
    table={};status={r:True for r in regions};widths=[]
    for r in intervals:
        key=(r['region'],r['cell'],r['gene'])
        if key in table or type(r['lower']) is not int or type(r['upper']) is not int or not 0<=r['lower']<=r['upper']:
            raise ValueError('duplicate or invalid count interval')
        table[key]=r;covered=r['lower']<=truth[key]<=r['upper']
        status[r['region']]=status.get(r['region'],True) and covered;widths.append(r['upper']-r['lower'])
    if not set(truth).issubset(table):raise ValueError('interval table omits nonzero true count entries')
    n=len(status);hits=sum(status.values())
    return {'RNA_records':len(records),'regions':n,'misassigned_RNA':sum(errors),
            'pooled_AUROC':misassignment_auc(scores,errors),
            'simultaneous_region_coverage':hits,'fraction':hits/n if n else None,
            'wilson95':wilson(hits,n) if n else None,
            'mean_interval_width':sum(widths)/len(widths) if widths else None,
            'scope':'Supplied aligned truth only; pooled AUROC has dependent within-region transcripts. Region confidence interval requires independent region units.'}


def fwer_gate(upper_pvalues, *, family_coverage_error, familywise_error):
    """Conditional statistical bridge, NOT an expression-to-p-value generator.

    Requires predeclared hypotheses, valid true-data null p-values, and simultaneous
    uncertainty-set coverage >=1-alpha. Each supplied upper p must bound the null
    p-value over the ENTIRE family. Then Bonferroni at (beta-alpha)/m controls
    unconditional FWER <= beta by a union bound, without independence assumptions.
    """
    alpha,beta=F(family_coverage_error),F(familywise_error)
    if not 0<=alpha<beta<1 or not upper_pvalues:raise ValueError('need 0<=coverage error<FWER<1 and nonempty hypotheses')
    pp=list(map(F,upper_pvalues))
    if any(p<0 or p>1 for p in pp):raise ValueError('p-value upper bounds must lie in [0,1]')
    cutoff=(beta-alpha)/len(pp)
    return {'threshold':str(cutoff),'reject':[p<=cutoff for p in pp],
            'guarantee':'FWER <= beta only under the documented coverage and valid-p-value assumptions; this function does not establish those premises'}
