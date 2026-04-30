"""Check p41270 value format (how multiple codes are stored in one field)."""
import subprocess, dxpy, os, csv

dispensed = dxpy.find_one_data_object(typename="Dataset", name="app*.dataset",
                                       folder="/", name_mode="glob")
dataset_id = dispensed["id"]

# Extract just eid + p41270 + p41280_a0 (to see format)
out_path = "/tmp/hes_format_test.csv"
r = subprocess.run(
    ["dx", "extract_dataset", dataset_id,
     "--fields", "participant.eid,participant.p41270,participant.p41280_a0",
     "-o", out_path, "--delimiter", ","],
    capture_output=True, text=True, timeout=200
)
print(f"rc={r.returncode}")
if r.returncode != 0:
    print("STDERR:", r.stderr[:500])
    raise SystemExit(1)

print(f"File size: {os.path.getsize(out_path)/1e6:.1f} MB")

# Show first 10 non-null rows
with open(out_path) as f:
    reader = csv.reader(f)
    header = next(reader)
    print(f"\nColumns: {header}")
    print("\nFirst 10 rows with non-null p41270:")
    count = 0
    for row in reader:
        if row[1] and row[1] not in ("", "[]", "null"):  # p41270 not empty
            print(f"  eid={row[0]}  p41270={row[1][:120]}  p41280_a0={row[2]}")
            count += 1
            if count >= 10:
                break
