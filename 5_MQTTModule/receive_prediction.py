import paho.mqtt.client as mqtt

# Callback when the client connects to the broker
def on_connect(client, userdata, flags, rc):
    print(f"Connected to MQTT broker with result code: {rc}")
    client.subscribe("roomTemp_predictions")

# Callback when a message is received from the subscribed topic
def on_message(client, userdata, msg):
    payload = msg.payload.decode('utf-8')
    print(f"Received: {payload}")

if __name__ == "__main__":
    # Create an MQTT client
    client = mqtt.Client()

    # Set the callback functions
    client.on_connect = on_connect
    client.on_message = on_message

    # Connect to the broker
    client.connect("mqtt.eclipseprojects.io", 1883, 60)

    # Start the MQTT loop to handle incoming messages
    client.loop_start()

    try:
        # Keep the script running until interrupted
        while True:
            pass

    except KeyboardInterrupt:
        # Disconnect the client when the script is interrupted
        client.disconnect()
        client.loop_stop()
