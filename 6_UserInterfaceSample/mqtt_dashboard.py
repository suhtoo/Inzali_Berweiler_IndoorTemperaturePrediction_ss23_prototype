import streamlit as st
import pandas as pd
import plotly.express as px

# Function to read data from the CSV file and create the dashboard
def create_dashboard():
    st.set_page_config(page_title="Indoor Temperature Prediction",
                   page_icon="🌡",
                   layout="wide")

    st.title("🌡 Indoor Temperature Predictions Dashboard")
    st.markdown("##")

    # Read the data from the CSV file
    data_df = pd.read_csv("received_data.csv")

    # Display the raw data table
    st.subheader("Raw Incoming Data")
    st.dataframe(data_df)
    st.text("Window Status : 1=closed, 2=tilted, 3=fully opened")

    # Create a line plot for room temperature over time
    fig = px.line(data_df, x='timestamp', y=['weatherTemp','roomTemp_predictions'], title='Room Temperature Over Time')
    st.plotly_chart(fig)

if __name__ == "__main__":
    create_dashboard()
