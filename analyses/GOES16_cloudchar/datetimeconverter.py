# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 12:37:45 2022

@author: Dean

Collection of functions to convert dates, times, and datetimes to various 
units.  
"""



import datetime



def yyyymmdd_to_ndays(date_str, with_year=False):
    """
    Convert date from a string 'yyyymmdd' to days since Jan1, 2020, +1 - i.e. Jan 1
    would not be 0, it would be 1. Default is to return result as an int; if 
    'with_year' is set to True, returns a string in the format 'yyyy###' where 
    '###' is the daynumber z-filled to 3 digits.
    
    Inputs:
        date_str: string '2020mmdd', m=month, d=day.
    """
    
    # Convert to datetime object, at T00:00:
    date_dt = datetime.datetime(
        year=int(date_str[0:4]), month=int(date_str[4:6]),
        day=int(date_str[6:8])
        )
    # Timedelta since Jan 1, 2020 - T00:00, +1 day:
    jan1_2020 = datetime.datetime(year=2020, month=1, day=1)
    daynum = (date_dt - jan1_2020).days + 1
    
    if with_year: return '2020' + str(daynum).zfill(3)
    return daynum
