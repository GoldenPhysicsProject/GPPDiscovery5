"""Write the leading n x n block of a dump (exact original strings) for s37. usage: subdump.py dump.json n out.json"""
import sys, json, hashlib
r = json.load(open(sys.argv[1])); n = int(sys.argv[2])
json.dump([x[:n] for x in r[:n]], open(sys.argv[3], 'w'))
print(hashlib.sha256(open(sys.argv[3], 'rb').read()).hexdigest())
