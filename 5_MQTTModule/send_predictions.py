import paho.mqtt.publish as publish
import time
import pandas as pd
from datetime import datetime, timedelta
import os

os.chdir('/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/5_MQTTModule')

# Function to publish roomTemp to the MQTT broker
def publish_roomTemp_to_mqtt(roomTemp_values):
    mqtt_broker = 'mqtt.eclipseprojects.io'
    mqtt_topic = "roomTemp_predictions"

    current_time = datetime.now()

    for i, roomTemp in enumerate(roomTemp_values):
        future_time = current_time + timedelta(minutes=5)
        current_time_str = current_time.strftime('%Y-%m-%d %H:%M:%S')
        payload = f"In 5 minutes, the predicted room temperature is {roomTemp:.2f}"
        publish.single(mqtt_topic, payload=payload, hostname=mqtt_broker)
        print(f"Sent: {payload}")
        time.sleep(30)  # To avoid flooding the MQTT broker

if __name__ == "__main__":
    # Load the predictions from the file (e.g., CSV)
    predictions_df = pd.read_csv("/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule/predictions.csv")
    roomTemp_predictions = predictions_df["roomTemp_predictions"].values

    # Send the predictions for 5 minutes into the future
    print("Sending roomTemp predictions for 5 minutes into the future...")
    publish_roomTemp_to_mqtt(roomTemp_predictions)
    print("Room temperature predictions sent for 5 minutes into the future!")
