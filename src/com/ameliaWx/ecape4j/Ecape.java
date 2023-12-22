package com.ameliaWx.ecape4j;

/**
 * 
 * @author Amelia Urquhart (https://github.com/a-urq)
 *
 */
public class Ecape {
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
		
		return updraftRadius;
	}
	
	/**
	 * 
	 * @param updraftRadius				Units: Meters
	 * @return <b>entrainmentRate:</b>	Units: m^-1
	 */
	public static double entrainmentRate(double updraftRadius) {
		double entrainmentRate = (2 * k2 * L_mix)/(Pr * updraftRadius * updraftRadius);
		
		return entrainmentRate;
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

	private static double specificHumidityFromDewpoint(double temperature, double dewpoint, double pressure) {
		double vaporPressure = vaporPressure(dewpoint);

		double specificHumidity = specificHumidity(pressure, vaporPressure, temperature);

		return specificHumidity;
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
	private static double specificHumidity(double pressure, double vaporPressure, double temperature) {
		double waterVaporDensity = absoluteHumidity(vaporPressure, temperature); // kg m^-3
		double airDensity = dryAirDensity(pressure - vaporPressure, temperature); // kg m^-3

		return waterVaporDensity / (waterVaporDensity + airDensity);
	}

	/** Units: J kg^-1 K^-1 */
	public static final double dryAirGasConstant = 287;
	/** Units: J kg^-1 K^-1 */
	public static final double waterVaporGasConstant = 461.5;
	/** Units: J kg^-1 */
	public static final double latentHeatOfVaporization = 2500000;

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
		double specificHumidity = specificHumidityFromDewpoint(temperature, dewpoint, pressure);

		return c_p * temperature + L_vr * specificHumidity + g * height;
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
