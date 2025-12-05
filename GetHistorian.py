import datetime as dt
import pymysql
import os
from pandas  import DataFrame

def GetHistorian(instrument, start_time, end_time=None,
                 hostname="sbcmirror.fnal.gov", port=3306,
                 user="reader", password_var="SQL_READER_PASSWORD",
                 database="SBCslowcontrol", table="DataStorage"):
    """
    Function that connects to a MySQL database, reads the value of instruments 
    between given timestamps, and returns the values in a pandas dataframe.

    Args:
      instrument (str): Name of the instrument to query.
      start_time (str or datetime.datetime): Starting timestamp.
      end_time (str or datetime.datetime, optional): Ending timestamp. 
        Defaults to current time if not provided.
      hostname (str, optional): Hostname of the MySQL server.
      port (int, optional): Port of the MySQL service.
      user (str, optional): MySQL user name.
      password_var (str, optional): Environment variable that saves the password.
        Set up by running `export SQL_READER_PASSWORD="samplepassword"`.
      database (str, optional): Name of the MySQL database.
      table (str, optional): Name of the table to query.
      
    Returns:
      pandas.DataFrame: DataFrame with columns ['Instrument', 'Time', 'Value'],
        containing the queried instrument data, sorted by time in ascending order.

    Raises:
      ValueError: If the password environment variable is not set.
      pymysql.Error: If the database connection or execution fails.
    """

    db = None
    cursor = None
    
    try: 
        # check parameters
        password = os.getenv(password_var)
        if password is None:
            raise ValueError(f"Password environment variable {password_var} is not declared.")
        if end_time is None: 
            end_time = dt.datetime.now()
    
        # open connection
        db = pymysql.connect(
                host=hostname, 
                user=user, 
                password=password,
                database=database, 
                port=port,
                connect_timeout=10)
        cursor = db.cursor()
    
        # prepare query
        query = f"""
        SELECT Instrument, Time, Value FROM {table} 
        WHERE Instrument = %s AND Time BETWEEN %s AND %s 
        ORDER BY Time ASC;
        """
        params = (instrument, start_time, end_time)
        
        cursor.execute(query, params)
        raw_data = cursor.fetchall()
        df = DataFrame(raw_data, columns=['Instrument', 'Time', 'Value'])

        return df
        
    # close connection
    finally:
        if cursor:
            cursor.close()
        if db:
            db.close()

if __name__ == "__main__":
    df = GetHistorian("PT1101", "2025-11-01 10:00", "2025-11-01 12:00")
    print(df.head())
