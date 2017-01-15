"""
Processing code for nextHIV
"""
from . import db, align, cluster, identifiers, resistance, phylogeny, sites, utils

def get_backend():
    return 'rethinkdb'
