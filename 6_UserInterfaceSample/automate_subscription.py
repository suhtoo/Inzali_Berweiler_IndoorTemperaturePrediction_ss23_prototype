import os
import time
import threading
import subprocess

# Function to run the MQTT subscriber and Streamlit dashboard
def run_subscriber():
    print("Running MQTT Subscriber...")
    subprocess.run(["python", "mqtt_subscriber.py"])

def run_dashboard():
    print("Running Streamlit Dashboard...")
    subprocess.run(["streamlit", "run", "mqtt_dashboard.py"])

if __name__ == "__main__":
    # Start the MQTT subscriber in a separate thread
    subscriber_thread = threading.Thread(target=run_subscriber)
    subscriber_thread.start()

    # Start the Streamlit dashboard
    run_dashboard()

