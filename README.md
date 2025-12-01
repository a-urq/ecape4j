# ecape4j
A simple Java library that computes ECAPE values and parcel paths. Code has been verified by checking it against the ECAPE_FUNCTIONS script written by Dr. Peters.

# Notice to Developers
Using this library on its own is quite boilerplate-heavy and may be unwieldy. The <a href="https://github.com/a-urq/weather-utils-java">WeatherUtils</a> ECAPE parcel method also uses this code, but will take care of the boilerplate on its own. I have decided to also upload this library on its own in case developers wish to use their own CAPE values, LFCs, ELs, storm motions, etc.

# Authors and Contributors
### **AUTHOR: Amelia R H Urquhart, University of Oklahoma** 
##### SPECIAL THANKS TO: Dr. John Peters, Pennsylvania State University

# Instructions for use
This library expects you to have already calculated a parcel origin state, storm motion vector, CAPE value, LFC height, and EL height. This is to allow for the user to supply their own values if they wish. To calculate these values, you may be interested in my other library <a href="https://github.com/a-urq/weather-utils-java">WeatherUtils</a>. 

The following code snippet demonstrates how Ecape4J may be used, making heavy use of WeatherUtils.

```java
// This code sets your input data
int numRecords = ...; // Number of records in your sounding data.

double[] pressure = new double[numRecords]; // Units: Pascals
double[] height = new double[numRecords]; // Units: Meters (both Above Sea Level and Above Ground Level will work)
double[] temperature = new double[numRecords]; // Units: Kelvins
double[] dewpoint = new double[numRecords]; // Units: Kelvins
double[] uWind = new double[numRecords]; // Units: m s^-1
double[] vWind = new double[numRecords]; // Units: m s^-1

// your code for entering the data into the arrays goes here
// remember to enter the lowest pressure (highest altitude) data in the first index and the highest pressure (lowest altitude) data in the last index.

// This code allows you to choose which CAPE type and storm motion vector you would like to use.
ParcelPath pathType = ParcelPath.MOST_UNSTABLE; // Valid options: SURFACE_BASED, MIXED_LAYER, MOST_UNSTABLE, INFLOW_LAYER
StormMotion smVector = StormMotion.BUNKERS_RIGHT; // Valid options: BUNKERS_RIGHT, BUNKERS_LEFT, MEAN_WIND

// This code computes the parcel origin state
ArrayList<RecordAtLevel> parcelPath = WeatherUtils.computeParcelPath(pressure, temperature, dewpoint, pathType,	false);
RecordAtLevel parcelOriginRaw = parcelPath.get(0);
double[] parcelOrigin = { parcelOriginRaw.pressure, parcelOriginRaw.height, parcelOriginRaw.temperature, parcelOriginRaw.dewpoint };

// This code computes CAPE, LFC, and EL based on the chosen CAPE type
double cape = WeatherUtils.computeCape(pressure, temperature, dewpoint, parcelPath);
double lfc = WeatherUtils.levelOfFreeConvection(pressure, temperature, dewpoint, parcelPath);
double el = WeatherUtils.equilibriumLevel(pressure, temperature, dewpoint, parcelPath);

// This code computes the storm motion vector
double[] stormMotion = null;
switch (smVector) {
case BUNKERS_LEFT:
	stormMotion = WeatherUtils.stormMotionBunkersIDLeftMoving(pressure, height, uWind, vWind);
	break;
case BUNKERS_RIGHT:
	stormMotion = WeatherUtils.stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);
	break;
case MEAN_WIND:
	stormMotion = WeatherUtils.stormMotionBunkersIDMeanWindComponent(pressure, height, uWind, vWind);
	break;
default:
	stormMotion = WeatherUtils.stormMotionBunkersIDRightMoving(pressure, height, uWind, vWind);
	break;
}

// This code specifies the inflow layer
double inflowBottom = 0; // Units: Meters
double inflowTop = 1000; // Units: Meters

boolean entrainmentSwitch = true;
boolean pseudoadiabaticSwitch = false;

double[][] ecapeParcelPathRaw = Ecape.ecapeParcel(pressure, height, temperature, dewpoint, uWind, vWind, entrainmentSwitch, pseudoadiabaticSwitch, parcelOrigin, stormMotion, inflowBottom, inflowTop, cape, lfc, el);
```

The value returned from `Ecape.ecapeParcel` is a two-dimensional array of doubles. Each element of the array is a list of pressures, heights, temperatures, water vapor mass fractions, or total water mass fractions. The following code snippet demonstrates how you might extract these lists. It is highly recommended that the user compute density temperature for plotting.

```java
double[] ecapeParcelPressure = ecapeParcelPathRaw[0]; // Units: Pascals
double[] ecapeParcelHeight = ecapeParcelPathRaw[1]; // Units: Meters
double[] ecapeParcelTemperature = ecapeParcelPathRaw[2]; // Units: Kelvins
double[] ecapeParcelQv = ecapeParcelPathRaw[3]; // Units: kg/kg
double[] ecapeParcelQt = ecapeParcelPathRaw[4]; // Units: kg/kg

double[] ecapeParcelDensityTemperature = Ecape.densityTemperature(ecapeParcelTemperature, ecapeParcelQv, ecapeParcelQt);
```
