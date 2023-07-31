import subprocess
import time
import os

#os.chdir('/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype')

def run_create_dataframe():
    subprocess.run(["python", "/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule/create_dataframe.py"])

def run_make_predictions():
    subprocess.run(["python", "/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule/make_predictions.py"])

def run_send_predictions():
    subprocess.run(["python", "/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/5_MQTTModule/send_predictions.py"])

if __name__ == "__main__":
    run_create_dataframe()
    time.sleep(5)  # Add a delay of 5 seconds if needed
    run_make_predictions()
    time.sleep(5)  # Add a delay of 5 seconds if needed
    run_send_predictions()
