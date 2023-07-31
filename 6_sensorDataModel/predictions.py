import pandas as pd
from pandas.plotting import scatter_matrix
import matplotlib.pyplot as plt
import sklearn as sk
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
import sensor_mqtt





def split_dataset():
    # read the dataset into a pandas dataframe
    df = pd.read_csv("C:/Users/Robin/Documents/ML/data/weatherDataAndRoomSensorData.csv")
    df = df.select_dtypes(exclude="object")
    
    integer_columns = df.select_dtypes(include='int64').columns
    df[integer_columns] = df[integer_columns].apply(pd.to_numeric)
    df.columns = df.columns.str.strip()
    df = df.drop('RTemperature (F)', axis=1)

    #print(df.columns)
    target = "RTemperature (C)" # only Temperature column
    X = df.drop(columns=target, axis=1)
    y = df[target]

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

    return X, y, X_train, X_test, y_train, y_test


def train_regression():
    x, y, X_train, X_test, y_train, y_test = split_dataset()
    model = LinearRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    print(y_pred)

def train_decision_tree():
    X, y, X_train, X_test, y_train, y_test = split_dataset()
    model = DecisionTreeRegressor()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    print(y_pred)

def train_random_forest():
    X, y, X_train, X_test, y_train, y_test = split_dataset()
    model = RandomForestRegressor()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    print(y_pred)

def fivefold_validation():
    X, y, X_train, X_test, y_train, y_test = split_dataset()
    model = RandomForestRegressor()
    model.fit(X_train, y_train)
    # Splitting the data into 5 folds
    kf = KFold(n_splits=5, shuffle=True)

    # Iterate over each fold
    for train_index, test_index in kf.split(X):

        # Split the data into training and testing sets
        X_train, X_test = X.iloc[train_index], X.iloc[test_index]
        y_train, y_test = y.iloc[train_index], y.iloc[test_index]

        # Train your model on the training set
        model.fit(X_train, y_train)

        # Test your model on the testing set
        accuracy = model.score(X_test, y_test)
        
        # Print the accuracy for this fold
        print("Accuracy:", accuracy)

n = 100

def current_data_test(n):
    
    df = sensor_mqtt.get_sensor_mqtt()
    print(df.shape)
    print("lalelu \n",df, df.shape)
    #df = 2023-07-27 22:35:38,70.00, 26.70, 80.06,16.8,91,1005,overcast clouds

    

    # data needs to be in the same format as test_data


    integer_columns = df.select_dtypes(include='int64').columns
    df[integer_columns] = df[integer_columns].apply(pd.to_numeric)
    df.columns = df.columns.str.strip()
    df = df.drop("RTemperature (C)", axis=1)
    df = df.drop("RTemperature (F)", axis=1)
    df = df.drop("Timestamp", axis=1)
    df = df.drop("Weather Conditions", axis=1)
    

    X, y, X_train, X_test, y_train, y_test = split_dataset()
    
    predictions = []
    for _ in range(n):
        model = DecisionTreeRegressor()
        model.fit(X_train, y_train)
        y_pred = model.predict(df)
        predictions.append(y_pred)

    pred_mean = sum(predictions) / len(predictions)
    print(pred_mean)
    return pred_mean


