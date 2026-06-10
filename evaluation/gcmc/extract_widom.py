import re
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import sys


root = Path(sys.argv[1])
pattern = re.compile(r"Average Henry coefficient:\s+([\d.e+\-]+)\s+\+/-\s+([\d.e+\-]+)")

records = []
for sub in tqdm(sorted(root.iterdir())):
    if not sub.is_dir():
        continue
    output_dir = sub / "output"
    if not output_dir.exists():
        continue
    for txt in output_dir.glob("*.txt"):
        text = txt.read_text(errors="ignore")
        m = pattern.search(text)
        if m:
            records.append({"name": sub.name, "KH": float(m.group(1)), "error": float(m.group(2))})
            break

df = pd.DataFrame(records, columns=["name", "KH", "error"])
out = Path(sys.argv[1]+"/widom_results_1.csv")
df.to_csv(out, index=False)
print(f"Saved {len(df)} records to {out}")
