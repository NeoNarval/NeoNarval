import numpy as np

""" 
This function will compare the two lists used as an input. It will return the proportion of values in the first list that find an equivalent with the given precision in the second list, and also a list with the matching values and their gap
"""
def data_matching(x,y,precision):
    
    matching_data = []
    
    for i in range(len(x)):
        gap_min = precision
        best_matching_value = None
        
        for j in range(len(y)):
                
                if ( abs(x[i]-y[j]) <= gap_min ):
                    
                    gap_min = abs(x[i]-y[j])
                    best_matching_value = y[j]
        data = [x[i],best_matching_value,gap_min]            
        if (best_matching_value != None ):
            matching_data.insert(i,data)
    
    matching_rate = float(len(matching_data))/float(len(x))
    print("Matched, Total",len(matching_data),len(x))
    return(matching_data,matching_rate)