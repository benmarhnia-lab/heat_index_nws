The .R file contains a function to calculate the heat index based on the National Weather Service modified Rothfusz regression, including the adjustments. I also calcualted the apparent temperature using the Kalkstein and Valimont 1986 equation. 

Input could be Fahrenheit/Celsius, dew point temperature/relative humidity; dew point temperature/relative humidity transformation is based on August-Roche-Magnus approximation.

Default output is a data.frame of four columns: heat index in Fahrenheit/Celsius and apparent temperature in Fahrenheit/Celsius, all replaced by air temperature when below 68 Fahrenheit.

Please contact Chen Chen (chc048@ucsd.edu) if have questions/suggestions.
