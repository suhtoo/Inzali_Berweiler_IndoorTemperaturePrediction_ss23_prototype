import pandas as pd
import time
import paho.mqtt.client as mqtt

# MQTT Broker details
broker_address = "mqtt.eclipseprojects.io"  # Replace with your MQTT broker address
broker_port = 1883  # Replace with your MQTT broker port
topic = "room_environment_data"  # The topic to which data will be published

def on_connect(client, userdata, flags, rc):
    if rc == 0:
        print("Connected to MQTT Broker!")
    else:
        print("Failed to connect, return code %d\n", rc)

def on_publish(client, userdata, mid):
    print("Data published successfully!")

def send_data_through_mqtt(dataframe):
    client = mqtt.Client()
    client.on_connect = on_connect
    client.on_publish = on_publish

    client.connect(broker_address, broker_port)
    client.loop_start()

    for _, row in dataframe.iterrows():
        values = row.values
        payload = ",".join(str(value) for value in values)
        print(f"Sending data: {payload}")
        client.publish(topic, payload)

        # Wait for 30 seconds
        time.sleep(30)

    client.loop_stop()
    client.disconnect()

if __name__ == "__main__":
    # Load the DataFrame from the CSV file
    predictions_df = pd.read_csv("/Users/hiz/Desktop/RWTH_3rdSem/CR_Protyping Project DSS - BIM/Project/Inzali_Berweiler_IndoorTemperaturePrediction_ss23_prototype/4_deploymentModule/predictions.csv")

    # Send the data through MQTT
    send_data_through_mqtt(predictions_df)
