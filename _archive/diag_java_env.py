"""Quick diagnostic: find Java location and test dx extract_dataset."""
import os, subprocess, glob, sys

print("=== ENVIRONMENT ===")
print(f"JAVA_HOME: {os.environ.get('JAVA_HOME', 'NOT SET')}")
print(f"PATH: {os.environ.get('PATH', '')}")
print(f"Python: {sys.executable}")

print("\n=== FIND JAVA BINARIES ===")
r = subprocess.run(["find", "/opt", "/usr", "/home", "-name", "java", "-type", "f",
                    "-maxdepth", "10"], capture_output=True, text=True, timeout=20)
print(r.stdout[:2000] or "(none found)")

print("\n=== JVM DIRECTORIES ===")
for pattern in ["/usr/lib/jvm/*/", "/opt/conda/pkgs/openjdk*/", "/opt/jdk*/", "/opt/java*/"]:
    found = glob.glob(pattern)
    if found:
        print(f"  {pattern} → {found}")

print("\n=== CONDA INFO ===")
r2 = subprocess.run(["conda", "info", "--json"], capture_output=True, text=True, timeout=10)
print(r2.stdout[:500] or r2.stderr[:200])

print("\n=== SOURCE CONDA + CHECK JAVA ===")
r3 = subprocess.run(
    'source /opt/conda/etc/profile.d/conda.sh 2>/dev/null; '
    'conda activate base 2>/dev/null; '
    'echo "JAVA_HOME=$JAVA_HOME"; '
    'which java 2>/dev/null || echo "java not in PATH"; '
    'java -version 2>&1 | head -2',
    shell=True, executable="/bin/bash", capture_output=True, text=True, timeout=30
)
print(r3.stdout[:500])
print(r3.stderr[:200])

print("\n=== TEST dx extract_dataset ===")
import dxpy
try:
    dispensed = dxpy.find_one_data_object(typename="Dataset", name="app*.dataset",
                                           folder="/", name_mode="glob")
    print(f"Dataset ID: {dispensed['id']}")
    r4 = subprocess.run(
        ["dx", "extract_dataset", dispensed["id"], "--list-fields"],
        capture_output=True, text=True, timeout=60
    )
    lines = r4.stdout.splitlines()
    print(f"  --list-fields returned {len(lines)} lines, rc={r4.returncode}")
    # Print first few lines and lines matching key fields
    print("  First 10 lines:")
    for l in lines[:10]:
        print(f"    {l}")
    print("  Lines matching 'eid':", [l for l in lines if "eid" in l.lower()][:3])
    print("  Lines matching 'p34':", [l for l in lines if " p34" in l or l.startswith("p34")][:3])
    print("  Lines matching 'p130000':", [l for l in lines if "p130000" in l][:3])
except Exception as e:
    print(f"dx extract_dataset test failed: {e}")

print("\n=== DONE ===")
