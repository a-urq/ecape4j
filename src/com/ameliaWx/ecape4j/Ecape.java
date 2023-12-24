package com.ameliaWx.ecape4j;

import java.util.ArrayList;

/**
 * 
 * @author Amelia Urquhart (https://github.com/a-urq)
 *
 */
public class Ecape {
//	public static void main(String[] args) {
////		double testPressure = 90000;
////		double testTemperature = 290;
////		double testDewpoint = 290;
////		double testHeight = 1000;
////
////		double specificHumidity = Ecape.specificHumidityFromDewpoint(testDewpoint, testPressure);
////		double specificHumidity2 = Ecape.specificHumidityFromDewpoint(testDewpoint, testPressure-10000);
////		
////		double resultDewpoint = Ecape.dewpointFromSpecificHumidity(specificHumidity, testPressure);
////		
////		System.out.println(specificHumidity);
////		System.out.println(specificHumidity2);
////		System.out.println(resultDewpoint);
//
////		double moistStaticEnergy = Ecape.moistStaticEnergy(testTemperature, testDewpoint, testHeight, testPressure);
//
////		double resultTemperature = Ecape.temperatureFromMoistStaticEnergy(moistStaticEnergy, testPressure, testHeight);
//
////		System.out.println(resultTemperature);
//		System.out.println(updraftRadius(4539, 2993, 10, 12000));
//		System.out.println(entrainmentRate(updraftRadius(1000, 800, 10, 10000)) + " m");
//	}

	/**
	 * Computes the value of Entrainment Cape as specified by John Peters et. al.
	 * 2023 (https://arxiv.org/pdf/2301.04712.pdf).
	 * 
	 * @param pressure           Units: Array of Pascals
	 * @param height             Units: Array of Meters
	 * @param temperature        Units: Array of Kelvins
	 * @param dewpoint           Units: Array of Kelvins
	 * @param uWind              Units: Array of m s^-1
	 * @param vWind              Units: Array of m s^-1
	 * @param stormMotion        Units: array of m s^-1, first entry: u-component,
	 *                           second entry: v-component
	 * @param inflowBottom       Units: Meters, lower bound of inflow layer
	 * @param inflowTop          Units: Meters, upper bound of inflow layer
	 * @param parcelOriginHeight Units: Meters
	 * @param cape               Units: J kg^-1
	 * @param lfc                Units: Meters
	 * @param el                 Units: Meters
	 * @return <b>entrainmentCape</b> Units: J kg^-1
	 */
	public static double entrainmentCape(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind, double[] stormMotion, double inflowBottom, double inflowTop,
			double parcelOriginHeight, double cape, double lfc, double el) {
		if (cape == 0) {
			return 0; // saves a lot of computing time
		}

		double psi = calcPsi(el);
		double vsr = calcVSR(height, uWind, vWind, stormMotion, inflowBottom, inflowTop);
		double[][] hHat = calcMseHat(pressure, height, temperature, dewpoint, parcelOriginHeight);
		double ncape = calcNcape(pressure, height, temperature, dewpoint, hHat, lfc, el);

		// terms separated out for readability and debuggability
		// making it all one line did NOT work out well
		double ecapeTerm1Numerator = -1 - psi - ((2 * psi) / (vsr * vsr)) * ncape;
		double ecapeTerm1And2Denominator = ((4 * psi) / (vsr * vsr)); // they are the same lol

		double radicandPart1 = 1 + psi + ((2 * psi) / (vsr * vsr)) * ncape;
		double radicandPart2 = ((8 * psi) / (vsr * vsr)) * (cape - psi * ncape);
		double ecapeTerm2Numerator = Math.sqrt(Math.pow(radicandPart1, 2) + radicandPart2);

		double ecape = ecapeTerm1Numerator / ecapeTerm1And2Denominator
				+ ecapeTerm2Numerator / ecapeTerm1And2Denominator;

//		double ecape = (-1 - psi - ((2*psi)/(vsr*vsr))*ncape)/((4*psi)/(vsr*vsr))
//				+ Math.pow(Math.pow(1 + psi + ((2*psi)/(vsr*vsr))*ncape + ((8*psi)/(vsr*vsr)) * (cape - psi*ncape), 2), 0.5)/
//				((4*psi)/(vsr*vsr));

		double ecapeA = vsr * vsr / 2.0 + ecape;

		if (ecapeA < 0) {
			ecapeA = 0;
		}

//		System.out.println("psi: " + psi);
//		System.out.println("vsr: " + vsr);
//		System.out.println("cape: " + cape);
//		System.out.println("ncape: " + ncape);
//		System.out.println("ecape: " + ecape);
//		System.out.println("ecape_a: " + ecapeA);

		return ecapeA;
	}

	/**
	 * 
	 * @param cape              Units: J kg^-1
	 * @param ecape             Units: J kg^-1
	 * @param vsr               Units: m s^-1
	 * @param stormColumnHeight Units: Meters
	 * @return <b>updraftRadius:</b> Units: Meters
	 */
	public static double updraftRadius(double cape, double ecape, double vsr, double stormColumnHeight) {
		double nondimE = ecape / cape;
		double nondimV = vsr / Math.sqrt(2 * ecape);

		double nondimR = Math
				.sqrt(((4 * sigma * sigma) / (alpha * alpha * Math.PI * Math.PI)) * ((nondimV * nondimV) / (nondimE)));

		double updraftRadius = nondimR * stormColumnHeight;

		return updraftRadius/2;
	}

	/**
	 * 
	 * @param updraftRadius Units: Meters
	 * @return <b>entrainmentRate:</b> Units: m^-1
	 */
	public static double entrainmentRate(double updraftRadius) {
		double entrainmentRate = (2 * k2 * L_mix) / (Pr * updraftRadius * updraftRadius);

//		entrainmentRate = 0.000046;
		
		return entrainmentRate;
	}

	private static final double ECAPE_PARCEL_DZ = 20.0; // Units: Meters
	private static final double DRY_ADIABATIC_LAPSE_RATE = 0.0098; // Units: K m^-1
	private static final double DEWPOINT_LAPSE_RATE = 0.0018; // Units: K m^-1

	/**
	 * 
	 * @param pressure
	 * @param height
	 * @param temperature
	 * @param dewpoint
	 * @param uWind
	 * @param vWind
	 * @param parcelOrigin { pressure [Pa], height [m], temperature [K], dewpoint
	 *                     [K] }
	 * @param stormMotion
	 * @param inflowBottom
	 * @param inflowTop
	 * @param cape
	 * @param lfc
	 * @param el
	 * @return <b>ecapeParcel:</b> 2D Array of Doubles <br>
	 *         Index 0 - Pressure [Pa] <br>
	 *         Index 1 - Height [m] <br>
	 *         Index 2 - Temperature [K] <br>
	 *         Index 3 - Dewpoint [K]
	 */
	public static double[][] ecapeParcel(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind, double[] parcelOrigin, double[] stormMotion, double inflowBottom,
			double inflowTop, double cape, double lfc, double el) {
		double ecape = entrainmentCape(pressure, height, temperature, dewpoint, uWind, vWind, stormMotion, inflowBottom,
				inflowTop, parcelOrigin[1], cape, lfc, el);
		double vsr = calcVSR(height, uWind, vWind, stormMotion, inflowBottom, inflowTop);

		double updraftRadius = updraftRadius(cape, ecape, vsr, el - parcelOrigin[1]);
		double entrainmentRate = entrainmentRate(updraftRadius);

		double[] moistStaticEnergy = new double[dewpoint.length];
		for (int i = 0; i < moistStaticEnergy.length; i++) {
			moistStaticEnergy[i] = moistStaticEnergy(temperature[i], dewpoint[i], height[i], pressure[i]);
		}

		double[] specificHumidity = new double[dewpoint.length];
		for (int i = 0; i < specificHumidity.length; i++) {
			double vaporPressure = Ecape.vaporPressure(dewpoint[i]);
			specificHumidity[i] = Ecape.specificHumidity(pressure[i], vaporPressure);
		}

		ArrayList<Double> pressureList = new ArrayList<>();
		ArrayList<Double> heightList = new ArrayList<>();
		ArrayList<Double> temperatureList = new ArrayList<>();
		ArrayList<Double> dewpointList = new ArrayList<>();
		ArrayList<Double> moistStaticEnergyList = new ArrayList<>();

		double parcelPressure = parcelOrigin[0];
		double parcelHeight = parcelOrigin[1];
		double parcelTemperature = parcelOrigin[2];
		double parcelDewpoint = parcelOrigin[3];

		double parcelMoistStaticEnergy = Ecape.moistStaticEnergy(parcelTemperature, parcelDewpoint, parcelHeight,
				parcelPressure);

		{
			pressureList.add(parcelPressure);
			heightList.add(parcelHeight);
			temperatureList.add(parcelTemperature);
			dewpointList.add(parcelDewpoint);
			moistStaticEnergyList.add(parcelMoistStaticEnergy);
		}

		while (parcelPressure > 10000) {
			parcelPressure = pressureAtHeight(parcelPressure, ECAPE_PARCEL_DZ, parcelTemperature);
			parcelHeight += ECAPE_PARCEL_DZ;

			if (parcelDewpoint < parcelTemperature) {
				parcelTemperature -= DRY_ADIABATIC_LAPSE_RATE * ECAPE_PARCEL_DZ;
				parcelDewpoint -= DEWPOINT_LAPSE_RATE * ECAPE_PARCEL_DZ;

				double parcelSpecificHumidity = specificHumidityFromDewpoint(parcelDewpoint, parcelPressure);
				double envSpecificHumidity = revLinearInterp(height, specificHumidity, parcelHeight);

				parcelSpecificHumidity += -entrainmentRate * (parcelSpecificHumidity - envSpecificHumidity)
						* ECAPE_PARCEL_DZ;

				parcelDewpoint = dewpointFromSpecificHumidity(parcelSpecificHumidity, parcelPressure);

				pressureList.add(parcelPressure);
				heightList.add(parcelHeight);
				temperatureList.add(parcelTemperature);
				dewpointList.add(parcelDewpoint);

				parcelMoistStaticEnergy = Ecape.moistStaticEnergy(parcelTemperature, parcelDewpoint, parcelHeight,
						parcelPressure);
				
//				System.out.printf("%8.2f Pa\tq=%8.6f\tq0=%8.6f\tdq=%8.6f\n", parcelPressure, parcelSpecificHumidity, envSpecificHumidity, -1000*entrainmentRate * (parcelSpecificHumidity - envSpecificHumidity)*ECAPE_PARCEL_DZ);
//				System.out.printf("%8.2f hPa\t%8.2f m\t%8.2f C\t%8.2f C\t%10.8f m\n", parcelPressure/100.0, parcelHeight, parcelTemperature - 273.15, parcelDewpoint - 273.15, updraftRadius);
			} else {
				double envMoistStaticEnergy = revLinearInterp(height, moistStaticEnergy, parcelHeight);

				parcelMoistStaticEnergy += -entrainmentRate * (parcelMoistStaticEnergy - envMoistStaticEnergy)
						* ECAPE_PARCEL_DZ;

				double parcelTemperatureFromMSE = Ecape.temperatureFromMoistStaticEnergy(parcelMoistStaticEnergy,
						parcelPressure, parcelHeight);
				
//				System.out.printf("%8.2f Pa\t%8.2f J kg^-1\t%8.2f J kg^-1\t%8.2f C\n", parcelPressure, parcelMoistStaticEnergy, envMoistStaticEnergy, parcelTemperatureFromMSE - 273.15);
//				System.out.printf("%8.2f hPa\t%8.2f m\t%8.2f C\t%8.2f C\t%10.8f m\t%10.8f m^-1\n", parcelPressure/100.0, parcelHeight, parcelTemperature - 273.15, parcelDewpoint - 273.15, updraftRadius, entrainmentRate);

				parcelTemperature = parcelTemperatureFromMSE;
				parcelDewpoint = parcelTemperatureFromMSE;
				
				pressureList.add(parcelPressure);
				heightList.add(parcelHeight);
				temperatureList.add(parcelTemperatureFromMSE);
				dewpointList.add(parcelTemperatureFromMSE);
			}
		}

		// ret[0]: pressure [Pa]
		// ret[1]: height [m]
		// ret[2]: temperature [K]
		// ret[3]: dewpoint [K]
		double[][] ret = new double[4][dewpointList.size()];

		for (int i = 0; i < ret[0].length; i++) {
			ret[0][i] = pressureList.get(i);
			ret[1][i] = heightList.get(i);
			ret[2][i] = temperatureList.get(i);
			ret[3][i] = dewpointList.get(i);
		}

		return ret;
	}

	private static final double c_p = 1005; // Units: J kg^-1
	private static final double L_vr = 2501000; // Units: J kg^-1
	private static final double g = 9.81; // Units: m s^-2
	private static final double sigma = 1.6; // Units: dimensionless
	private static final double alpha = 0.8; // Units: dimensionless
	private static final double k2 = 0.18; // Units: dimensionless
	private static final double Pr = 1.0 / 3.0; // Units: dimensionless
	private static final double L_mix = 120.0; // Units: Meters

	/**
	 * 
	 * @param equilibriumLevel Units: Meters
	 * @return <b>Psi:</b> For use in ECAPE computations
	 */
	private static double calcPsi(double equilibriumLevel) {
		double numerator = k2 * alpha * alpha * Math.PI * Math.PI * L_mix;
		double denominator = Pr * sigma * sigma * equilibriumLevel;

		return numerator / denominator;
	};

	/**
	 * 
	 * @param height
	 * @param uWind
	 * @param vWind
	 * @param stormMotion
	 * @param inflowBottom
	 * @param inflowTop
	 * @return <b>VSR:</b> For use in ECAPE computations
	 */
	private static double calcVSR(double[] height, double[] uWind, double[] vWind, double[] stormMotion,
			double inflowBottom, double inflowTop) {
		// DEFAULTS TO 0-1 KM INFLOW IF NO EIL IS FOUND
		if(inflowTop <= inflowBottom) {
			inflowBottom = 0;
			inflowTop = 1000;
		}
		 
		double uWindSum = 0.0;
		double vWindSum = 0.0;
		double weightSum = 0.0;

		for (int i = 0; i < uWind.length - 1; i++) {
			double heightUp = height[i];
			double heightDn = height[i + 1];

			if (heightUp > inflowTop && heightDn <= inflowTop) {
				double weight = (inflowTop - height[i + 1]);

				double uWindUp = revLinearInterp(height, uWind, heightUp);
				double vWindUp = revLinearInterp(height, vWind, heightUp);

				double uWindDn = uWind[i + 1];
				double vWindDn = vWind[i + 1];

				uWindSum += weight * (uWindUp + uWindDn) / 2.0;
				vWindSum += weight * (vWindUp + vWindDn) / 2.0;
				weightSum += weight;
			} else if (heightDn > inflowTop && heightUp <= inflowBottom) {
				double weight = (height[i] - height[i + 1]);

				double uWindUp = uWind[i];
				double vWindUp = vWind[i];

				double uWindDn = uWind[i + 1];
				double vWindDn = vWind[i + 1];

				uWindSum += weight * (uWindUp + uWindDn) / 2.0;
				vWindSum += weight * (vWindUp + vWindDn) / 2.0;
				weightSum += weight;
			} else if (heightUp > inflowBottom && heightDn <= inflowBottom) {
				double weight = (height[i] - inflowBottom);

				double uWindUp = uWind[i];
				double vWindUp = vWind[i];

				double uWindDn = revLinearInterp(height, uWind, heightDn);
				double vWindDn = revLinearInterp(height, vWind, heightDn);

				uWindSum += weight * (uWindUp + uWindDn) / 2.0;
				vWindSum += weight * (vWindUp + vWindDn) / 2.0;
				weightSum += weight;
			}
		}

		double uWindAvg = uWindSum / weightSum;
		double vWindAvg = vWindSum / weightSum;

		return Math.hypot(uWindAvg - stormMotion[0], vWindAvg - stormMotion[1]);
	}

	private static final double MSEHAT_DZ = 20; // Units: Meters

	/**
	 * mseHat[0] = function input (Units: Meters) mseHat[1] = function output
	 * (Units: J kg^-1)
	 * 
	 * @param height
	 * @param temperature
	 * @param dewpoint
	 * @return
	 */
	private static double[][] calcMseHat(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double parcelOriginHeight) {
		double heightTop = height[0];

		int numInputs = (int) ((heightTop - parcelOriginHeight) / MSEHAT_DZ);

		double[] moistStaticEnergy = new double[dewpoint.length];
		for (int i = 0; i < moistStaticEnergy.length; i++) {
			moistStaticEnergy[i] = moistStaticEnergy(temperature[i], dewpoint[i], height[i], pressure[i]);
		}

		double[][] ret = new double[2][numInputs];

		double mseRunningAvg = 0.0;
		for (int i = 0; i < numInputs; i++) {
			double mseAtZ = revLinearInterp(height, moistStaticEnergy, parcelOriginHeight + i * MSEHAT_DZ);

			mseRunningAvg = (i) / (i + 1.0) * mseRunningAvg + 1.0 / (i + 1.0) * mseAtZ;

			ret[0][i] = parcelOriginHeight + i * MSEHAT_DZ;
			ret[1][i] = mseRunningAvg;
		}

		return ret;
	}

	private static final double NCAPE_DZ = 20; // Units: Meters

	/**
	 * Not to be confused with normalized CAPE, this is just for the ECAPE
	 * calculation
	 * 
	 * @param pressure
	 * @param height
	 * @param temperature
	 * @param dewpoint
	 * @param hHat
	 * @param lfc
	 * @param el
	 * @return
	 */
	private static double calcNcape(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[][] hHat, double lfc, double el) {
		double[] hStar = new double[dewpoint.length];
		for (int i = 0; i < hStar.length; i++) {
			hStar[i] = moistStaticEnergy(temperature[i], temperature[i], height[i], pressure[i]);
		}

		double ncape = 0.0;

		for (double z = lfc; z < el; z += NCAPE_DZ) {
			double hHatAtZ = linearInterp(hHat[0], hHat[1], z);
			double hStarAtZ = revLinearInterp(height, hStar, z);

			double temperatureAtZ = revLinearInterp(height, temperature, z);

			double integrand = -g / (c_p * temperatureAtZ) * (hHatAtZ - hStarAtZ);

			ncape += integrand * NCAPE_DZ;
		}

		return ncape;
	}

	private static double specificHumidityFromDewpoint(double dewpoint, double pressure) {
		double vaporPressure = vaporPressure(dewpoint);

		double specificHumidity = specificHumidity(pressure, vaporPressure);

		return specificHumidity;
	}

	private static double dewpointFromSpecificHumidity(double specificHumidity, double pressure) {
		double vaporPressure = vaporPressureFromSpecificHumidity(pressure, specificHumidity);

		double dewpoint = dewpointFromVaporPressure(vaporPressure);

		return dewpoint;
	}

	private static double vaporPressureFromSpecificHumidity(double pressure, double specificHumidity) {
		double numerator = specificHumidity * pressure;
		double denominatorTerm = (1 / waterVaporGasConstant + specificHumidity / dryAirGasConstant
				- specificHumidity / waterVaporGasConstant);

		double vaporPressure = numerator / (dryAirGasConstant * denominatorTerm);

		return vaporPressure;
	}

	private static double dewpointFromVaporPressure(double vaporPressure) {
		double e0 = 611; // Pascals
		double t0 = 273.15; // Kelvins

		double dewpointReciprocal = 1 / t0
				- waterVaporGasConstant / latentHeatOfVaporization * Math.log(vaporPressure / e0);

		return 1 / dewpointReciprocal;
	}

	/**
	 * Computes the specific humidity using total pressure, vapor pressure, and
	 * temperature.
	 * 
	 * @param pressure      Units: Pascals
	 * @param vaporPressure Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>specificHumidity</b> Units: Fraction
	 */
	private static double specificHumidity(double pressure, double vaporPressure) {
		double waterVaporDensity = absoluteHumidity(vaporPressure, 280); // kg m^-3
		double airDensity = dryAirDensity(pressure - vaporPressure, 280); // kg m^-3

		return waterVaporDensity / (waterVaporDensity + airDensity);
	}

	/** Units: J kg^-1 K^-1 */
	private static final double dryAirGasConstant = 287;
	/** Units: J kg^-1 K^-1 */
	private static final double waterVaporGasConstant = 461.5;
	/** Units: J kg^-1 */
	private static final double latentHeatOfVaporization = 2500000;

	/**
	 * Computes the partial density of water vapor, also called absolute humidity.
	 * 
	 * @param vaporPressure Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>absoluteHumidity</b> Units: kg m^-3
	 */
	private static double absoluteHumidity(double vaporPressure, double temperature) {
		double waterVaporDensity = vaporPressure / (waterVaporGasConstant * temperature);

		return waterVaporDensity;
	}

	/**
	 * Computes the partial density of dry air.
	 * 
	 * @param dryAirPressure Units: Pascals
	 * @param temperature    Units: Kelvins
	 * @return <b>dryAirDensity</b> Units: kg m^-3
	 */
	private static double dryAirDensity(double dryAirPressure, double temperature) {
		double dryAirDensity = dryAirPressure / (dryAirGasConstant * temperature);

		return dryAirDensity;
	}

	/**
	 * Computes the partial pressure of water vapor in air with the dewpoint given.
	 * 
	 * @param dewpoint Units: Kelvins
	 * @return <b>vaporPressure</b> Units: Pascals
	 */
	private static double vaporPressure(double dewpoint) {
		double e0 = 611; // Pascals
		double t0 = 273.15; // Kelvins

		return e0 * Math.exp(latentHeatOfVaporization / waterVaporGasConstant * (1 / t0 - 1 / dewpoint));
	}

	private static double moistStaticEnergy(double temperature, double dewpoint, double height, double pressure) {
		double specificHumidity = specificHumidityFromDewpoint(dewpoint, pressure);

		return c_p * temperature + L_vr * specificHumidity + g * height;
	}

	/**
	 * Assumes parcel is SATURATED, decodes MSE into temperature from ECAPE parcel
	 */
	private static double temperatureFromMoistStaticEnergy(double mse, double p, double z) {
		double moistNonstaticEnergy = mse - z * g;
		double guessT = moistNonstaticEnergy / c_p;
		double guessMSE = Ecape.moistStaticEnergy(guessT, guessT, z, p);

		final double guessChangeCoef = 0.2; // lmao this is entirely arbitrary
		while (Math.abs(mse - guessMSE) > 1) {
			double guessDiff = mse - guessMSE;

			guessT -= -guessDiff / c_p * guessChangeCoef;
			guessMSE = Ecape.moistStaticEnergy(guessT, guessT, z, p);
		}

		return guessT;
	}

//	/**
//	 * Computes the mixing ratio using total pressure, vapor pressure, and
//	 * temperature.
//	 * 
//	 * Borrowed from WeatherUtils.
//	 * 
//	 * @param pressure Units: Pascals
//	 * @param dewpoint Units: Kelvins
//	 * @return <b>mixingRatio</b> Units: Fraction (kg kg^-1)
//	 */
//	private static double mixingRatio(double pressure, double dewpoint) {
//		double vaporPressure = vaporPressure(dewpoint);
//		double mixingRatio = 0.62197 * (vaporPressure) / (pressure - vaporPressure);
//
//		return mixingRatio;
//	}

	/** Units: J K^-1 mol^-1 */
	private static final double molarGasConstant = 8.314;
	/** Units: kg mol^-1 */
	private static final double avgMolarMass = 0.029;

	/**
	 * Computes pressure at a given height above sea level. Uses temperature to
	 * figure out scale height.
	 * 
	 * @param seaLevelPres        Units: Pascals
	 * @param heightAboveSeaLevel Units: Meters
	 * @param temperature         Units: Kelvins
	 * @return <b>pressureAtHeight</b> Units: Pascals
	 */
	private static double pressureAtHeight(double seaLevelPres, double heightAboveSeaLevel, double temperature) {
		double scaleHeight = (molarGasConstant * temperature) / (avgMolarMass * g); // Meters

		return seaLevelPres * Math.exp(-heightAboveSeaLevel / scaleHeight);
	}

	// inputArr assumed to already be sorted and increasing
	private static double linearInterp(double[] inputArr, double[] outputArr, double input) {
		if (input < inputArr[0]) {
			return outputArr[0];
		} else if (input >= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for (int i = 0; i < inputArr.length - 1; i++) {
				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];

				if (input == input1) {
					return outputArr[i];
				} else if (input < input2) {
					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];

					double weight1 = (input2 - input) / (input2 - input1);
					double weight2 = (input - input1) / (input2 - input1);

					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}

			return -1024.0;
		}
	}

	// inputArr assumed to already be sorted and decreases
	private static double revLinearInterp(double[] inputArr, double[] outputArr, double input) {
		if (input > inputArr[0]) {
			return outputArr[0];
		} else if (input <= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for (int i = 0; i < inputArr.length - 1; i++) {
				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];

				if (input == input2) {
					return outputArr[i];
				} else if (input > input2) {
					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];

					double weight1 = (input2 - input) / (input2 - input1);
					double weight2 = (input - input1) / (input2 - input1);

					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}

			return -1024.0;
		}
	}
}
