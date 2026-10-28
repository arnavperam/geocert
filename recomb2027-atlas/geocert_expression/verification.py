"""Bind a penalized flow proof to the requested model and molecular score.

This independently checks flow and physical witnesses, but recompiles using the
production compiler. It is not a proof-assistant verification of that compiler.
"""
from collections import Counter
from fractions import Fraction as F
from .budget import check_scores
from pinned_compiler import graph_for
from check_flow import check_certificate
from mixed import verify_mixed

def verify_penalized(engine,scores,answer):
    scores=check_scores(engine.problem,scores)
    if answer.get('query_key')!=engine.score_key(scores):raise ValueError('wrong query identity')
    mu=F(answer['multiplier']);sign=answer['sign']
    if mu<0 or sign not in (-1,1):raise ValueError('invalid objective')
    cert=answer.get('certificate')
    if cert is None:raise ValueError('flow proof missing')
    graph=graph_for(engine.compiled,[[-sign*v for v in row] for row in scores],mu)
    if (cert['n'],cert['s'],cert['t'],cert['scale'],F(cert['constant'])) != (graph['n'],graph['s'],graph['t'],graph['scale'],graph['constant']):raise ValueError('graph metadata mismatch')
    if Counter(tuple(row[:3]) for row in cert['arcs'])!=Counter(graph['edges']):raise ValueError('capacity graph mismatch')
    check_certificate(cert)
    feasible=cert['value']<graph['big']
    if feasible!=answer['feasible']:raise ValueError('wrong feasibility flag')
    if not feasible:return True
    side=set(cert['source_side']);owners=engine.problem.decode([int(i in side) for i in range(engine.problem.n*engine.problem.m)])
    if owners!=answer['owners']:raise ValueError('owners do not match certified cut')
    score=sum(row[o] for row,o in zip(scores,owners));area=F(answer['area'])
    if score!=answer['score'] or mu*area-sign*score!=F(answer['objective']):raise ValueError('physical objective mismatch')
    if F(answer['objective'])!=F(cert['value'],cert['scale'])+F(cert['constant']):raise ValueError('proof objective mismatch')
    verify_mixed(engine.problem,owners,answer['boundaries'],area)
    return True
