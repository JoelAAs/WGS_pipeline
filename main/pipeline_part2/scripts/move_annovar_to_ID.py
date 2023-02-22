import sys
import re

rs_pattern = re.compile("snp138=(rs[0-9]+);")
id_pattern = re.compile("\t(\\.)\t")
for line in sys.stdin:
	if line[0] == "#":
		print(line, end="")
	else:
		pos_match = rs_pattern.search(line)
		if pos_match:
			id_idx = id_pattern.search(line).span(1)
			found_rs = pos_match.group(1)
			print(line[:id_idx[0]], end="")
			print(found_rs, end="")
			print(line[id_idx[1]:], end="")
		else:
			print(line, end="")
