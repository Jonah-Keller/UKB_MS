"""Find actual column names for HES fields p41270 and p41280."""
import subprocess, dxpy

dispensed = dxpy.find_one_data_object(typename="Dataset", name="app*.dataset",
                                       folder="/", name_mode="glob")
dataset_id = dispensed["id"]

r = subprocess.run(["dx", "extract_dataset", dataset_id, "--list-fields"],
                   capture_output=True, text=True, timeout=120)
lines = r.stdout.splitlines()
print(f"Total fields: {len(lines)}, rc={r.returncode}")

# Print all lines mentioning p41270 or p41280
hes_lines = [l for l in lines if "p41270" in l or "p41280" in l]
print(f"\nHES field lines ({len(hes_lines)} found):")
for l in hes_lines[:30]:
    print(" ", l)

# Also check p34 (year_birth) naming
p34_lines = [l for l in lines if "p34\t" in l or "p34_" in l or l.startswith("participant.p34")]
print(f"\np34 lines:")
for l in p34_lines[:5]:
    print(" ", l)
