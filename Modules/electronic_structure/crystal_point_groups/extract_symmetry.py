import json
from lxml import html
import requests

point_groups = {}

for n in range(1,33):
	page = requests.get('http://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-point_wp-list?num='+str(n))
	tree = html.fromstring(page.content)
	name = tree[1][2].text.replace("Wyckoff Positions of the 3D Crystallographic Point Group ","")
	for t in tree[1][2]:
		if t.tag == "sub":
			name += "_"+t.text
		if t.tail is not None:
			name += t.tail
	axes = None
	if " [Unique axes b]" in name:
		name = name.replace(" [Unique axes b]","")
		axes = "Unique axes b"
	if " [Hexagonal axes]" in name:
		name = name.replace(" [Hexagonal axes]", "")
		axes = "Hexagonal axes"
	Shoenflies, Hermann = name.split("(")
	Hermann = Hermann.replace(")", "").strip()

	ops = tree[1][4][0][4][3].text
	ops = ops.replace(")(","),(")
	ops = ops.replace("(","{").replace(")","}")
	ops = ops.replace(",",", ")
	ops = "{" + ops +"}"
	
	point_groups[n] = {"ops": ops, "Hermann": Hermann, "Shoenflies": Shoenflies, "axes": axes}

json.dump(point_groups, open("crystal_point_groups.json", "w"), indent=2)
