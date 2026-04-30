import subprocess
import time
import sys

JOB_ID = "job-J7P959QJpjp7159VJj3PyZxP"

def check_job():
    while True:
        res = subprocess.run(
            ["/Users/jonahkeller/ELab/CADASIL_Proteome_ML_Keller_2024_Rebuttal/.venv/bin/dx", "describe", JOB_ID, "--json"],
            capture_output=True, text=True
        )
        if res.returncode != 0:
            print(f"Error checking job: {res.stderr}")
            time.sleep(15)
            continue
            
        import json
        data = json.loads(res.stdout)
        state = data.get("state")
        
        print(f"Job {JOB_ID} state: {state}")
        
        if state == "done":
            print("\nJob completed successfully! Integration script starting...")
            # Run the integration script
            download = subprocess.run([sys.executable, "analysis/00_rap_extraction/04_download_and_integrate.py"])
            if download.returncode == 0:
                print("Data is successfully on disk!")
            break
        elif state in ["failed", "terminated", "partially_failed"]:
            print(f"\nJob failed with state: {state}")
            break
            
        time.sleep(30)

if __name__ == "__main__":
    check_job()
