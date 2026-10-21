"""Grok S67-resumed: extract T from an s66 mori_route pickle dump into the JSON row format read by s36_cholres.py."""
import sys, pickle, json
d = pickle.load(open(sys.argv[1], 'rb')); json.dump(d['T'], open(sys.argv[2], 'w'))
