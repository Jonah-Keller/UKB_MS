import subprocess, time, sys, json
JOB_ID = "job-J7P959QJpjp7159VJj3PyZxP"

while True:
    res = subprocess.run(
        ["/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx", "describe", JOB_ID, "--json"],
        capture_output=True, text=True
    )
    if res.returncode == 0:
        state = json.loads(res.stdout).get("state")
        if state not in ["runnable", "running", "waiting_on_input"]:
            print(f"Final state: {state}")
            break
    time.sleep(30)
