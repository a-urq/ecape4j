package com.ameliaWx.ecape4j;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * 
 * @author Amelia Urquhart (https://github.com/a-urq)
 *
 */

// GET METPY FORMULAS FOR SPECIFIC HUMIDITY/DEWPOINT CONVERSION
// THOSE HERE ARE KNOWN TO BE INCORRECT IN UPPER-ATMOSPHERIC CONDITIONS
public class Ecape {
	// Constant Declarations
	private static final double R_d = 287.04; // J kg^-1 K^-1
	private static final double R_v = 461.5; // J kg^-1 K^-1
	private static final double c_pd = 1005; // Units: J kg^-1
	private static final double c_pv = 1870; // Units: J kg^-1
	private static final double c_pl = 4190; // Units: J kg^-1
	private static final double c_pi = 2106; // Units: J kg^-1
	private static final double L_v_trip = 2501000; // Units: J kg^-1
	private static final double L_i_trip = 333000; // Units: J kg^-1
	private static final double T_trip = 273.15; // Units: J kg^-1

	private static final double g = 9.81; // Units: m s^-2
	private static final double sigma = 1.1; // Units: dimensionless
	private static final double alpha = 0.8; // Units: dimensionless
	private static final double k2 = 0.18; // Units: dimensionless
	private static final double Pr = 1.0 / 3.0; // Units: dimensionless
	private static final double L_mix = 120.0; // Units: Meters

	public static void main(String[] args) {
		double testPressure = 90000;
		double testTemperature = 290;
		double testDewpoint = 290;
		double testHeight = 1000;

		double specificHumidity = Ecape.specificHumidityFromDewpoint(testDewpoint, testPressure);
		double specificHumidity2 = Ecape.specificHumidityFromDewpoint(305, 100000);

		double resultDewpoint = Ecape.dewpointFromSpecificHumidity(specificHumidity, testPressure);
		
		System.out.println("saturatedAdiabaticLapseRate: " + saturatedAdiabaticLapseRate(238, 0.005, 35000, 233, 0.0001, 0, 0, -1024));

//		System.out.println(specificHumidity);
//		System.out.println(specificHumidity2);
//		System.out.println(resultDewpoint);

//		double moistStaticEnergy = Ecape.moistStaticEnergy(testTemperature, testDewpoint, testHeight, testPressure);

//		double resultTemperature = Ecape.temperatureFromMoistStaticEnergy(moistStaticEnergy, testPressure, testHeight);

//		System.out.println(updraftRadius(3558, 872, 3.44, 15004));
//		System.out.println(entrainmentRate(updraftRadius(4233.88, 3911.55, 18.17, 12500)) + " m");
	}

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

//		System.out.println("pressure: " + Arrays.toString(pressure));
//		System.out.println("height: " + Arrays.toString(height));
//		System.out.println("temperature: " + Arrays.toString(temperature));
//		System.out.println("dewpoint: " + Arrays.toString(dewpoint));
//		System.out.println("hHat[0]: " + Arrays.toString(hHat[0]));
//		System.out.println("hHat[1]: " + Arrays.toString(hHat[1]));

		// terms separated out for readability and debuggability
		// making it all one line did NOT work out well
		double ecapeTerm1Numerator = -1 - psi - ((2 * psi) / (vsr * vsr)) * ncape;
		double ecapeTerm1And2Denominator = ((4 * psi) / (vsr * vsr)); // they are the same lol

		double radicandPart1 = 1 + psi + ((2 * psi) / (vsr * vsr)) * ncape;
		double radicandPart2 = ((8 * psi) / (vsr * vsr)) * (cape - psi * ncape);
		double ecapeTerm2Numerator = Math.sqrt(Math.pow(radicandPart1, 2) + radicandPart2);

		double ecape = ecapeTerm1Numerator / ecapeTerm1And2Denominator
				+ ecapeTerm2Numerator / ecapeTerm1And2Denominator;
		double ecapeA = vsr * vsr / 2.0 + ecape;

		if (ecapeA < 0) {
			ecapeA = 0;
		}

		return ecapeA;
	}

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
	 * @return {<b>entrainmentCape</b> Units: J kg^-1, <b>ncape</b> Units: J kg^-1}
	 */
	public static double[] ecapeAndNcape(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind, double[] stormMotion, double inflowBottom, double inflowTop,
			double parcelOriginHeight, double cape, double lfc, double el) {
		if (cape == 0) {
			return new double[] { 0, 0 }; // saves a lot of computing time
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
		double ecapeA = vsr * vsr / 2.0 + ecape;

		if (ecapeA < 0) {
			ecapeA = 0;
		}

		return new double[] { ecapeA, ncape };
	}

	/**
	 * 
	 * @param cape              Units: J kg^-1
	 * @param ecape             Units: J kg^-1
	 * @param ncape             Units: J kg^-1
	 * @param vsr               Units: m s^-1
	 * @param stormColumnHeight Units: Meters
	 * @return <b>entrainmentRate:</b> Units: m^-1
	 */
	public static double entrainmentRate(double cape, double ecape, double ncape, double vsr,
			double stormColumnHeight) {
		double eaTilde = ecape / cape;
		double nTilde = ncape / cape;
		double vsrTilde = vsr / Math.sqrt(2 * cape);

		double eTilde = eaTilde - (vsrTilde * vsrTilde);

//		System.out.println("cape: " + cape);
//		System.out.println("ecape: " + ecape);
//		System.out.println("ncape: " + ncape);
//		System.out.println(eaTilde);
//		System.out.println(nTilde);
//		System.out.println(vsrTilde);

		double entrainmentRate = ((2 * (1 - eTilde)) / (eTilde + nTilde)) / (stormColumnHeight);

		return entrainmentRate;
	}

	/**
	 * 
	 * @param entrainmentRate Units: m^-1
	 * @return <b>updraftRadius:</b> Units: Meters
	 */
	public static double updraftRadius(double entrainmentRate) {
		double updraftRadius = Math.sqrt(2 * k2 * L_mix / (Pr * entrainmentRate));

		return updraftRadius;
	}

	/**
	 * 
	 * @param cape              Units: J kg^-1
	 * @param ecape             Units: J kg^-1
	 * @param vsr               Units: m s^-1
	 * @param stormColumnHeight Units: Meters
	 * @return <b>updraftRadius:</b> Units: Meters
	 */
	public static double updraftRadiusLegacy(double cape, double ecape, double vsr, double stormColumnHeight) {
		double nondimE = ecape / cape;
		double nondimV = vsr / Math.sqrt(2 * ecape);

		double nondimR = Math
				.sqrt(((4 * sigma * sigma) / (alpha * alpha * Math.PI * Math.PI)) * ((nondimV * nondimV) / (nondimE)));

		double updraftRadius = nondimR * stormColumnHeight;

		return updraftRadius / 2;
	}

	/**
	 * 
	 * @param updraftRadius Units: Meters
	 * @return <b>entrainmentRate:</b> Units: m^-1
	 */
	public static double entrainmentRateLegacy(double updraftRadius) {
		double entrainmentRate = (2 * k2 * L_mix) / (Pr * updraftRadius * updraftRadius);

//		entrainmentRate = 0.000046;

		return entrainmentRate;
	}

	private static final double ECAPE_PARCEL_DZ = 20.0; // Units: Meters
	private static final double DRY_ADIABATIC_LAPSE_RATE = 0.0098; // Units: K m^-1
	private static final double DEWPOINT_LAPSE_RATE = 0.0018; // Units: K m^-1

	/**
	 * Despite the name, it computes all four types of parcels from Dr Peters's
	 * formulas, including
	 * 
	 * - Pseudoadiabatic Non-entraining (traditional method, but does include more nuances
	 * - Pseudoadiabatic Entraining
	 * - Irreversible Adiabatic Non-entraining
	 * - Irreversible Adiabatic Entraining (recommended method)
	 * 
	 * @param pressure
	 * @param height
	 * @param temperature
	 * @param dewpoint
	 * @param uWind
	 * @param vWind
	 * @param entrainmentSwitch
	 * @param pseudoadiabaticSwitch
	 * @param parcelOrigin          { pressure [Pa], height [m], temperature [K],
	 *                              dewpoint [K] }
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
	 *         Index 3 - Water Vapor Mass Fraction [kg/kg] <br>
	 *         Index 4 - Total Water Mass Fraction (Vapor + Condensate) [kg/kg]
	 */
	public static double[][] ecapeParcel(double[] pressure, double[] height, double[] temperature, double[] dewpoint,
			double[] uWind, double[] vWind, boolean entrainmentSwitch, boolean pseudoadiabaticSwitch,
			double[] parcelOrigin, double[] stormMotion, double inflowBottom, double inflowTop, double cape, double lfc,
			double el) {
		double[] specificHumidity = new double[dewpoint.length];
		
		for(int i = 0; i < specificHumidity.length; i++) {
			specificHumidity[i] = specificHumidityFromDewpoint(dewpoint[i], pressure[i]); // REMEMBER THAT THIS NEEDS TO BE FIXED
		}
		
		if((cape == -1024 || lfc == -1024 || el == -1024) && entrainmentSwitch == true) {
			double[][] undilutedParcelProfile = ecapeParcel(pressure, height, temperature, dewpoint, uWind, vWind, false, pseudoadiabaticSwitch, parcelOrigin, stormMotion, inflowBottom, inflowTop, cape, lfc, el);

			double[] undilProfZ = undilutedParcelProfile[1];
			double[] undilProfT = undilutedParcelProfile[2];
			double[] undilProfQv = undilutedParcelProfile[3];
			double[] undilProfQt = undilutedParcelProfile[4];
			
			double[] capeLfcEl = computeCapeLfcEl(undilProfZ, undilProfT, undilProfQv, undilProfQt, height, temperature, specificHumidity);
			
//			System.out.println("undilProfZ:");
//			for(int i = 0; i < undilProfZ.length; i++) {
//				System.out.print(undilProfZ[i] + ", ");
//				
//				if(i % 10 == 0) System.out.print("<" + undilutedParcelProfile[0][i] + ">\n");
//			}
			
			if(cape == -1024) {
				cape = capeLfcEl[0];
			}
			
			if(lfc == -1024) {
				lfc = capeLfcEl[1];
			}
			
			if(el == -1024) {
				el = capeLfcEl[2];
			}
		}
		
		double entrainmentRate = 0.0;
		
		if (entrainmentSwitch) {
//			System.out.println("pressure: " + Arrays.toString(pressure));
//			System.out.println("stormMotion: " + Arrays.toString(stormMotion));
//			System.out.println("inflowBottom: " + inflowBottom);
//			System.out.println("inflowTop: " + inflowTop);
//			System.out.println("parcelOriginHeight: " + parcelOrigin[1]);
//			System.out.println("cape: " + cape);
//			System.out.println("lfc: " + lfc);
//			System.out.println("el: " + el);
			
			double[] ecapeNcape = ecapeAndNcape(pressure, height, temperature, dewpoint,
					uWind, vWind, stormMotion, inflowBottom, inflowTop,
					parcelOrigin[1], cape, lfc, el);
			double vsr = calcVSR(height, uWind, vWind, stormMotion, inflowBottom, inflowTop);

//			System.out.println(Arrays.toString(ecapeNcape));

			entrainmentRate = entrainmentRate(cape, ecapeNcape[0], ecapeNcape[1], vsr, el - parcelOrigin[1]);

//			System.out.println("Entrainment rate - CAPE: " + cape);
//			System.out.println("Entrainment rate - ECAPE: " + ecapeNcape[0]);
//			System.out.println("Entrainment rate - NCAPE: " + ecapeNcape[1]);
//			System.out.println("Entrainment rate - VSR: " + vsr);
//			System.out.println("Entrainment rate - COLUMN: " + (el - parcelOrigin[1]));
//			System.out.println("Entrainment rate: " + (entrainmentRate));
//			System.out.println("Updraft radius: " + updraftRadius(entrainmentRate));
		}
		
		double parcelPressure = parcelOrigin[0];
		double parcelHeight = logInterp(pressure, height, parcelPressure);
		double parcelTemperature = parcelOrigin[2];
		double parcelQv = specificHumidityFromDewpoint(parcelOrigin[3], parcelPressure);
		double parcelQt = parcelQv;
		
//		System.out.println("parcel origin pressure: " + parcelPressure);
//		System.out.println("parcel origin height: " + parcelHeight);
//		System.out.println("parcel origin temperature: " + parcelTemperature);
//		System.out.println("parcel origin qv: " + parcelQv);
//		System.out.println("parcel origin qt: " + parcelQt);
//		System.out.println("entrainment rate: " + entrainmentRate);

		ArrayList<Double> pressureList = new ArrayList<>();
		ArrayList<Double> heightList = new ArrayList<>();
		ArrayList<Double> temperatureList = new ArrayList<>();
		ArrayList<Double> qvList = new ArrayList<>();
		ArrayList<Double> qtList = new ArrayList<>();
		
		pressureList.add(parcelPressure);
		heightList.add(parcelHeight);
		temperatureList.add(parcelTemperature);
		qvList.add(parcelQv);
		qtList.add(parcelQt);
		
		double prate = 1 / ECAPE_PARCEL_DZ;
		
		if (!pseudoadiabaticSwitch) {
			prate = 0;
		}
		
		double dqtDz = 0;
		
		while(parcelPressure > pressure[0] && parcelPressure >= 10000) {
			double envTemperature = revLinearInterp(height, temperature, parcelHeight + parcelOrigin[1]);
			
			double parcelSaturationQv = (1 - parcelQt) * rSat(parcelTemperature, parcelPressure, 1);
			
			if(parcelSaturationQv > parcelQv) {
				parcelPressure = pressureAtHeight(parcelPressure, ECAPE_PARCEL_DZ, envTemperature);
				parcelHeight += ECAPE_PARCEL_DZ;

				envTemperature = revLinearInterp(height, temperature, parcelHeight);
				double envQv = revLinearInterp(height, specificHumidity, parcelHeight);
				
				double dTdz = unsaturatedAdiabaticLapseRate(parcelTemperature, parcelQv, envTemperature, envQv, entrainmentRate);
				double dqvDz = -entrainmentRate * (parcelQv - envQv);
				
				parcelTemperature += dTdz * ECAPE_PARCEL_DZ;
				parcelQv += dqvDz * ECAPE_PARCEL_DZ;
				parcelQt = parcelQv;
			} else {
				parcelPressure = pressureAtHeight(parcelPressure, ECAPE_PARCEL_DZ, envTemperature);
				parcelHeight += ECAPE_PARCEL_DZ;

				envTemperature = revLinearInterp(height, temperature, parcelHeight);
				double envQv = revLinearInterp(height, specificHumidity, parcelHeight);
				
				double dTdz = -1024;
				
				if(pseudoadiabaticSwitch) {
					dTdz = saturatedAdiabaticLapseRate(parcelTemperature, parcelQt, parcelPressure, envTemperature, envQv, entrainmentRate, prate, dqtDz);
				} else {
					dTdz = saturatedAdiabaticLapseRate(parcelTemperature, parcelQt, parcelPressure, envTemperature, envQv, entrainmentRate, prate, -1024.0);
				}
				
				double newParcelQv = (1 - parcelQt) * rSat(parcelTemperature, parcelPressure, 1);
				
				if(pseudoadiabaticSwitch) {
					dqtDz = (newParcelQv - parcelQv) / ECAPE_PARCEL_DZ;
				} else {
					dqtDz = -entrainmentRate * (parcelQt - envQv) - prate * (parcelQt - parcelQv);
				}
				
				parcelTemperature += dTdz * ECAPE_PARCEL_DZ;
				parcelQv = newParcelQv;
				
//				System.out.println("dTdz: " + dTdz);
//				System.out.println("z: " + parcelHeight);
//				System.out.println("parcelTemperature-sat: " + parcelTemperature);
				
				if(Double.isNaN(parcelTemperature)) break;
				
				if (pseudoadiabaticSwitch) {
					parcelQt = parcelQv;
				} else {
					dqtDz = -entrainmentRate * (parcelQt - envQv) - prate * (parcelQt - parcelQv);
					parcelQt += dqtDz * ECAPE_PARCEL_DZ;
				}
				
				if(parcelQt < parcelQv) {
					parcelQv = parcelQt;
				}
			}

//			System.out.println("adding to list");
			pressureList.add(parcelPressure);
			heightList.add(parcelHeight);
			temperatureList.add(parcelTemperature);
			qvList.add(parcelQv);
			qtList.add(parcelQt);
//			System.out.println("temperatureList.size(): " + temperatureList.size());
		}

		// ret[0]: pressure [Pa]
		// ret[1]: height [m]
		// ret[2]: temperature [K]
		// ret[3]: water vapor mass fraction [kg/kg]
		// ret[4]: total water mass fraction (vapor + condensate) [kg/kg]
		double[][] ret = new double[5][temperatureList.size()];
		
//		System.out.println("ret size: " + ret[0].length);

		for (int i = 0; i < ret[0].length; i++) {
			ret[0][i] = pressureList.get(i);
			ret[1][i] = heightList.get(i);
			ret[2][i] = temperatureList.get(i);
			ret[3][i] = qvList.get(i);
			ret[4][i] = qtList.get(i);
		}

		return ret;
	}
	
	private static double[] computeCapeLfcEl(double[] zParcel, double[] tParcel, double[] qvParcel, double[] qtParcel, double[] zEnv, double[] tEnv, double[] qvEnv) {
		double[] tRhoParcel = densityTemperature(tParcel, qvParcel, qtParcel);
		double[] tRhoEnv = densityTemperature(tEnv, qvEnv, qvEnv);
		
		double integratedPositiveBuoyancy = 0.0;
		double lfc = -1024;
		double el = -1024;
		
//		double[] mseEnv = moistStaticEnergy(zEnv, tEnv, qvEnv);

//		System.out.println("mseEnv: " + Arrays.toString(mseEnv));
		
//		int mseMinIdx = indexOfMinimum(mseEnv);
//		double zMinMse = zEnv[mseMinIdx];
//		System.out.println("mseMinIdx: " + Arrays.toString(mseEnv));
		
		for(int i = zParcel.length - 1; i >= 1; i--) {
			double z0 = zParcel[i];
			double dz = zParcel[i] - zParcel[i - 1];
			
			double tRho0 = revLinearInterp(zEnv, tRhoEnv, z0);
			double tRho = tRhoParcel[i];
			
			double buoyancy = g * (tRho - tRho0) / tRho0;
			
//			System.out.println(tRho + "\t" + tRho0 + "\t" + buoyancy + "\t" + z0 + "\t" + zMinMse);
			
			if (buoyancy > 0 && el == -1024) {
				el = z0;
			}
			
			if (buoyancy > 0 && lfc == -1024) {
				integratedPositiveBuoyancy += buoyancy * dz;
			}
			
			if (buoyancy < 0 && el != -1024 && z0 < 5000) {
				lfc = z0;
				break;
			}
		}
		
		if(lfc == -1024) {
			lfc = zParcel[0];
		}
		
		return new double[] {integratedPositiveBuoyancy, lfc, el};
	}
	
	private static int indexOfMinimum(double[] arr) {
		double min = Double.MAX_VALUE;
		int minIdx = -1;
		
		for(int i = 0; i < arr.length; i++) {
			if(arr[i] < min) {
				min = arr[i];
				minIdx = i;
			}
		}
		
		return minIdx;
	}

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
	public static double[][] ecapeParcelLegacy(double[] pressure, double[] height, double[] temperature,
			double[] dewpoint, double[] uWind, double[] vWind, double[] parcelOrigin, double[] stormMotion,
			double inflowBottom, double inflowTop, double cape, double lfc, double el) {
		long ecapeParcelStart = System.currentTimeMillis();

		double ecape = entrainmentCape(pressure, height, temperature, dewpoint, uWind, vWind, stormMotion, inflowBottom,
				inflowTop, parcelOrigin[1], cape, lfc, el);
		double vsr = calcVSR(height, uWind, vWind, stormMotion, inflowBottom, inflowTop);

		double updraftRadius = updraftRadiusLegacy(cape, ecape, vsr, el - parcelOrigin[1]);
		double entrainmentRate = entrainmentRateLegacy(updraftRadius);

//		System.out.println("updraft radius: " + updraftRadius);
//		System.out.println("entrainment rate: " + entrainmentRate);

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

		while (parcelPressure > pressure[0] && ecape > 0) {
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

				if (entrainmentRate > 1 / ECAPE_PARCEL_DZ) {
					entrainmentRate = 0.95 * 1 / ECAPE_PARCEL_DZ; // keeps numerical errors from going unstable. it's
																	// not like it's gonna matter much if the
																	// entrainment rate is already bonkers high
				}

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

		long ecapeParcelEnd = System.currentTimeMillis();

		System.out.println("ecape parcel compute time: " + (ecapeParcelEnd - ecapeParcelStart) + " ms");

		return ret;
	}

	// Equation 19 in Peters et. al. 2022
	// (https://journals.ametsoc.org/view/journals/atsc/79/3/JAS-D-21-0118.1.xml)
	public static double unsaturatedAdiabaticLapseRate(double temperatureParcel, double qvParcel, double temperatureEnv,
			double qvEnv, double entrainmentRate) {
		double temperatureEntrainment = -entrainmentRate * (temperatureParcel - temperatureEnv);

		double densityTemperatureParcel = densityTemperature(temperatureParcel, qvParcel, qvParcel);
		double densityTemperatureEnv = densityTemperature(temperatureEnv, qvEnv, qvEnv);

		double buoyancy = g * (densityTemperatureParcel - densityTemperatureEnv) / densityTemperatureEnv;

		double c_pmv = (1 - qvParcel) * c_pd + qvParcel * c_pv;

		double term1 = -g / c_pd;
		double term2 = 1 + (buoyancy / g);
		double term3 = c_pmv / c_pd;

		double dTdz = term1 * (term2 / term3) + temperatureEntrainment;

		return dTdz;
	}

	// Equation 24 in Peters et. al. 2022
	// (https://journals.ametsoc.org/view/journals/atsc/79/3/JAS-D-21-0118.1.xml)
	public static double saturatedAdiabaticLapseRate(double temperatureParcel, double qtParcel, double pressureParcel,
			double temperatureEnv, double qvEnv, double entrainmentRate, double prate, double qtEntrainment) {
		return saturatedAdiabaticLapseRate(temperatureParcel, qtParcel, pressureParcel, temperatureEnv, qvEnv,
				entrainmentRate, prate, qtEntrainment, 273.15, 253.15);
	}

	// Equation 24 in Peters et. al. 2022
	// (https://journals.ametsoc.org/view/journals/atsc/79/3/JAS-D-21-0118.1.xml)
	public static double saturatedAdiabaticLapseRate(double temperatureParcel, double qtParcel, double pressureParcel,
			double temperatureEnv, double qvEnv, double entrainmentRate, double prate, double qtEntrainment,
			double warmestMixedPhaseTemp, double coldestMixedPhaseTemp) {
//		System.out.println(temperatureParcel + "\t" + qtParcel + "\t" + pressureParcel + "\t" + temperatureEnv + "\t" + qvEnv + "\t" + entrainmentRate + "\t" + prate + "\t" + qtEntrainment);
		
		double omega = iceFraction(temperatureParcel, warmestMixedPhaseTemp, coldestMixedPhaseTemp);
		double dOmega = iceFractionDeriv(temperatureParcel, warmestMixedPhaseTemp, coldestMixedPhaseTemp);

		double qVsl = (1 - qtParcel) * rSat(temperatureParcel, pressureParcel, 0);
		double qVsi = (1 - qtParcel) * rSat(temperatureParcel, pressureParcel, 2);

		double qvParcel = (1 - omega) * qVsl + omega * qVsi;

		double temperatureEntrainment = -entrainmentRate * (temperatureParcel - temperatureEnv);
		double qvEntrainment = -entrainmentRate * (qvParcel - qvEnv);

		if (qtEntrainment == -1024.0) {
			qtEntrainment = -entrainmentRate * (qtParcel - qvEnv) - prate * (qtParcel - qvParcel);
		}

		double qCondensate = qtParcel - qvParcel;

		double qlParcel = qCondensate * (1 - omega);
		double qiParcel = qCondensate * omega;

		double c_pm = (1 - qtParcel) * c_pd + qvParcel * c_pv + qlParcel * c_pl + qiParcel * c_pi;

		double densityTemperatureParcel = densityTemperature(temperatureParcel, qvParcel, qtParcel);
		double densityTemperatureEnv = densityTemperature(temperatureEnv, qvEnv, qvEnv);

		double buoyancy = g * (densityTemperatureParcel - densityTemperatureEnv) / densityTemperatureEnv;

		double L_v = L_v_trip + (temperatureParcel - T_trip) * (c_pv - c_pl);
		double L_i = L_i_trip + (temperatureParcel - T_trip) * (c_pl - c_pi);

		double L_s = L_v + omega * L_i;

		double Q_vsl = qVsl / (PHI - PHI * qtParcel + qvParcel);
		double Q_vsi = qVsi / (PHI - PHI * qtParcel + qvParcel);

		double Q_M = (1 - omega) * (qVsl) / (1 - Q_vsl) + omega * (qVsi) / (1 - Q_vsi);
		double L_M = (1 - omega) * L_v * (qVsl) / (1 - Q_vsl) + omega * (L_v + L_i) * (qVsi) / (1 - Q_vsi);
		double R_m0 = (1 - qvEnv) * R_d + qvEnv * R_v;

		double term1 = buoyancy;
		double term2 = g;
		double term3 = ((L_s * Q_M) / (R_m0 * temperatureEnv)) * g;
		double term4 = (c_pm - L_i * (qtParcel - qvParcel) * dOmega) * temperatureEntrainment;
		double term5 = L_s * (qvEntrainment + qvParcel / (1 - qtParcel) * qtEntrainment);

		double term6 = c_pm;
		double term7 = (L_i * (qtParcel - qvParcel) - L_s * (qVsi - qVsl)) * dOmega;
		double term8 = (L_s * L_M) / (R_v * temperatureParcel * temperatureParcel);
		
//		System.out.println(L_s);
//		System.out.println(qvEntrainment);
//		System.out.println(qvParcel);
//		System.out.println(qtParcel);
//		System.out.println(qtEntrainment);
//		
//		System.out.println(term1);
//		System.out.println(term2);
//		System.out.println(term3);
//		System.out.println(term4);
//		System.out.println(term5);
//		System.out.println(term6);
//		System.out.println(term7);
//		System.out.println(term8);
		
		return -(term1 + term2 + term3 - term4 - term5) / (term6 - term7 + term8);
	}

	private static final double PHI = R_d / R_v;

	public static double densityTemperature(double temperature, double qv, double qt) {
		double densityTemperature = temperature * (1 - qt + qv / PHI);

		return densityTemperature;
	}

	public static double[] densityTemperature(double[] temperature, double[] qv, double[] qt) {
		double[] densityTemperature = new double[temperature.length];
		
		for(int i = 0; i < qv.length; i++) {
			densityTemperature[i] = densityTemperature(temperature[i], qv[i], qt[i]);
		}
		
		return densityTemperature;
	}

	private static double iceFraction(double temperature, double warmestMixedPhaseTemp, double coldestMixedPhaseTemp) {
		if (temperature >= warmestMixedPhaseTemp) {
			return 0;
		} else if (temperature <= coldestMixedPhaseTemp) {
			return 1;
		} else {
			return (1 / (coldestMixedPhaseTemp - warmestMixedPhaseTemp)) * (temperature - warmestMixedPhaseTemp);
		}
	}

	private static double iceFractionDeriv(double temperature, double warmestMixedPhaseTemp,
			double coldestMixedPhaseTemp) {
		if (temperature >= warmestMixedPhaseTemp) {
			return 0;
		} else if (temperature <= coldestMixedPhaseTemp) {
			return 0;
		} else {
			return (1 / (coldestMixedPhaseTemp - warmestMixedPhaseTemp));
		}
	}

	private static double rSat(double temperature, double pressure, int iceFlag) {
		return rSat(temperature, pressure, iceFlag, 273.15, 253.15);
	}

	private static final double vaporPresRef = 611.2;

	private static double rSat(double temperature, double pressure, int iceFlag, double warmestMixedPhaseTemp,
			double coldestMixedPhaseTemp) {
		switch (iceFlag) {
		case 0:
			double term1 = (c_pv - c_pl) / R_v;
			double term2 = (L_v_trip - T_trip * (c_pv - c_pl)) / R_v;

			double esi = Math.exp((temperature - T_trip) * term2 / (temperature * T_trip)) * vaporPresRef
					* Math.pow(temperature / T_trip, term1);

			double qSat = PHI * esi / (pressure - esi);

			return qSat;
		case 1:
			double omega = iceFraction(temperature, warmestMixedPhaseTemp, coldestMixedPhaseTemp);

			double qSat_l = rSat(temperature, pressure, 0, warmestMixedPhaseTemp, coldestMixedPhaseTemp);
			double qSat_i = rSat(temperature, pressure, 2, warmestMixedPhaseTemp, coldestMixedPhaseTemp);

			qSat = (1 - omega) * qSat_l + omega * qSat_i;

			return qSat;
		case 2:
			term1 = (c_pv - c_pi) / R_v;
			term2 = (L_v_trip - T_trip * (c_pv - c_pi)) / R_v;

			esi = Math.exp((temperature - T_trip) * term2 / (temperature * T_trip)) * vaporPresRef
					* Math.pow(temperature / T_trip, term1);

			qSat = PHI * esi / (pressure - esi);

			return qSat;
		default:
			return rSat(temperature, pressure, iceFlag, warmestMixedPhaseTemp, coldestMixedPhaseTemp);
		}
	}

	/**
	 * 
	 * @param equilibriumLevel Units: Meters
	 * @return <b>Psi:</b> For use in ECAPE computations
	 */
	private static double calcPsi(double equilibriumLevel) {
		double numerator = k2 * alpha * alpha * Math.PI * Math.PI * L_mix;
		double denominator = 4 * Pr * sigma * sigma * equilibriumLevel;

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
		if (inflowTop <= inflowBottom) {
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
//		System.out.println("height: " + Arrays.toString(height));

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
			
//			System.out.printf("%5.0f m\t%f.0 J/kg\t%f.0 J/kg\n", ret[0][i], ret[1][i], mseAtZ);
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
		
//		System.out.println("NCAPE-LFC: " + lfc + " m");
//		System.out.println("NCAPE-EL: " + el + " m");

		for (double z = lfc; z < el; z += NCAPE_DZ) {
			double hHatAtZ = linearInterp(hHat[0], hHat[1], z);
			double hStarAtZ = revLinearInterp(height, hStar, z);

			double temperatureAtZ = revLinearInterp(height, temperature, z);

			double integrand = -g / (c_pd * temperatureAtZ) * (hHatAtZ - hStarAtZ);

			ncape += integrand * NCAPE_DZ;
			
//			if(z % 500 == 0) System.out.printf("%5.0f m\t%f5.1 K\t%f.0 J/kg\t%f.0 J/kg\t%f.0 J/kg\t%f.0 J/kg\n", z, temperatureAtZ, hHatAtZ, hStarAtZ, integrand, ncape);
		}
		
//		System.out.println("Final NCAPE: " + ncape + " J/kg");

		return ncape;
	}

	private static double specificHumidityFromDewpoint(double dewpoint, double pressure) {
		double vaporPressure = vaporPresRef * Math.exp(17.67*(dewpoint - 273.15)/(dewpoint - 29.65));
		
		double mixingRatio = PHI * vaporPressure / (pressure - vaporPressure);

		double specificHumidity = mixingRatio / (1 + mixingRatio);

		return specificHumidity;
	}

	private static double dewpointFromSpecificHumidity(double specificHumidity, double pressure) {
		double mixingRatio = specificHumidity / (1 - specificHumidity);
		
		double vaporPressure = mixingRatio * pressure / (PHI + mixingRatio);
		
		double dewpoint = (273.15 - Math.log(vaporPressure/vaporPresRef)*29.65/17.67)/(1 - Math.log(vaporPressure/vaporPresRef)/17.67);

		return dewpoint;
	}

	private static double vaporPressureFromSpecificHumidity(double pressure, double specificHumidity) {
		double numerator = specificHumidity * pressure;
		double denominatorTerm = (1 / R_v + specificHumidity / R_d - specificHumidity / R_v);

		double vaporPressure = numerator / (R_d * denominatorTerm);

		return vaporPressure;
	}

	private static double dewpointFromVaporPressure(double vaporPressure) {
		double e0 = 611; // Pascals
		double t0 = 273.15; // Kelvins

		double dewpointReciprocal = 1 / t0 - R_v / L_v_trip * Math.log(vaporPressure / e0);

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

//		System.out.println("d_wv: " + waterVaporDensity);
//		System.out.println("d_da: " + airDensity);
//		System.out.println("q: " + waterVaporDensity / (waterVaporDensity + airDensity));

		return waterVaporDensity / (waterVaporDensity + airDensity);
	}

	/**
	 * Computes the partial density of water vapor, also called absolute humidity.
	 * 
	 * @param vaporPressure Units: Pascals
	 * @param temperature   Units: Kelvins
	 * @return <b>absoluteHumidity</b> Units: kg m^-3
	 */
	private static double absoluteHumidity(double vaporPressure, double temperature) {
		double waterVaporDensity = vaporPressure / (R_v * temperature);

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
		double dryAirDensity = dryAirPressure / (R_d * temperature);

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

		return e0 * Math.exp(L_v_trip / R_v * (1 / t0 - 1 / dewpoint));
	}

	private static double moistStaticEnergy(double temperature, double dewpoint, double height, double pressure) {
		double specificHumidity = specificHumidityFromDewpoint(dewpoint, pressure);

		return moistStaticEnergy(temperature, specificHumidity, height);
	}

	private static double moistStaticEnergy(double temperature, double specificHumidity, double height) {
		return c_pd * temperature + L_v_trip * specificHumidity + g * height;
	}

	public static double[] moistStaticEnergy(double[] temperature, double[] specificHumidity, double[] height) {
		double[] moistStaticEnergy = new double[specificHumidity.length];
		
		System.out.println();
		
		for(int i = 0; i < specificHumidity.length; i++) {
			moistStaticEnergy[i] = moistStaticEnergy(temperature[i], specificHumidity[i], height[i]);
		}
		
		return moistStaticEnergy;
	}

	/**
	 * Assumes parcel is SATURATED, decodes MSE into temperature from ECAPE parcel
	 */
	private static double temperatureFromMoistStaticEnergy(double mse, double p, double z) {
		double moistNonstaticEnergy = mse - z * g;
		double guessT = moistNonstaticEnergy / c_pd;
		double guessMSE = Ecape.moistStaticEnergy(guessT, guessT, z, p);

		final double guessChangeCoef = 0.2; // lmao this is entirely arbitrary
		while (Math.abs(mse - guessMSE) > 1) {
			double guessDiff = mse - guessMSE;

			guessT -= -guessDiff / c_pd * guessChangeCoef;
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

	// inputArr assumed to already be sorted and increasing
	private static double logInterp(double[] inputArr, double[] outputArr, double input) {
		if (input < inputArr[0]) {
			return outputArr[0];
		} else if (input >= inputArr[inputArr.length - 1]) {
			return outputArr[outputArr.length - 1];
		} else {
			for (int i = 0; i < inputArr.length - 1; i++) {
				if (i + 1 == outputArr.length) {
					return outputArr[outputArr.length - 1];
				}

				double input1 = inputArr[i];
				double input2 = inputArr[i + 1];

				if (input == input1) {
					return outputArr[i];
				} else if (input < input2) {
					double logInput1 = Math.log(input1);
					double logInput2 = Math.log(input2);
					double logInput = Math.log(input);

					double output1 = outputArr[i];
					double output2 = outputArr[i + 1];

					double weight1 = (logInput2 - logInput) / (logInput2 - logInput1);
					double weight2 = (logInput - logInput1) / (logInput2 - logInput1);

					return output1 * weight1 + output2 * weight2;
				} else {
					continue;
				}
			}

			return -1024.0;
		}
	}
}
