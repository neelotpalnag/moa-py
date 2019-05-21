package jmetal.init;

import java.util.HashMap;
import java.util.Properties;

import jmetal.core.Operator;

public class Params {
	
	static Properties properties=new Properties();
	static Operator selectionOperator;
	static Operator crossoverOperator;
	
	static double[][] hypervolumeTrueFront;
  static long startTime;
	
  public static void setStartTime(long time)
	{
  	startTime = time;
	}
	public static long getStartTime()
	{
		return startTime;
	}
	
	public static void setProperties(Properties p)
	{
		properties = p;
	}
	
	public static void setSeed(double seed)
	{
		properties.setProperty("Seed", String.valueOf(seed));
	}
	public static double getSeed()
	{
		return Double.parseDouble(properties.getProperty("Seed"));
	}
	
	public static void setSelectionOperator(Operator select)
	{
		selectionOperator = select;
	}
	public static Operator getSelectionOperator()
	{
		return selectionOperator;
	}
	
	public static void setCrossoverOperator(Operator crossover)
	{
		crossoverOperator = crossover;
	}
	public static Operator getCrossoverOperator()
	{
		return crossoverOperator;
	}
	
	public static int getHypervolumeSwitch()
	{
		return Integer.parseInt(properties.getProperty("Output.hypervolume"));
	}
	
	public static int getNumberOfObjectives()
	{
		return Integer.parseInt(properties.getProperty("Problem.numberOfObjectives"));
	}
	
	public static int getNumberOfVariables()
	{
		return Integer.parseInt(properties.getProperty("Problem.numberOfVariables"));
	}
	
	public static int getNumberOfFeatures()
	{
		return Integer.parseInt(properties.getProperty("Problem.numberOfFeatures"));
	}
	
	public static String getSimulatorSimIn()
	{
		return properties.getProperty("Simulator.simIn");
	}
	
	public static String getSimulatorSimOut()
	{
		return properties.getProperty("Simulator.simOut");
	}
	
	public static String getSimulatorCommand()
	{
		return properties.getProperty("Simulator.command");
	}
	
	public static double[] getVariablesLowerLimit()
	{
		int noOfVar = Integer.parseInt(properties.getProperty("Problem.numberOfVariables"));
		double [] limits = new double[noOfVar];
		String stringLimits = properties.getProperty("Variables.lowerLimit");
		String [] stringTableOfLimits = stringLimits.split(" ");
		int counter =0;
		
		for(int i=0; i< stringTableOfLimits.length; i++)
			if(stringTableOfLimits[i] != "" && stringTableOfLimits[i] != " " && stringTableOfLimits[i].length() > 0)
			{
				limits[counter] = Double.parseDouble(stringTableOfLimits[i]);
				counter++;
			}
		
		return limits;
	}
	
	public static double[] getVariablesUpperLimit()
	{
		int noOfVar = Integer.parseInt(properties.getProperty("Problem.numberOfVariables"));
		double [] limits = new double[noOfVar];
		String stringLimits = properties.getProperty("Variables.upperLimit");
		String [] stringTableOfLimits = stringLimits.split(" ");
		int counter =0;
		
		for(int i=0; i< stringTableOfLimits.length; i++)
		{
			if(stringTableOfLimits[i] != ""  && stringTableOfLimits[i] != " " && stringTableOfLimits[i].length() > 0)
			{
				limits[counter] = Double.parseDouble(stringTableOfLimits[i]);
				counter++;
			}
		}
		
		return limits;
	}

	public static void setDiscretizationStep(double[] ds) {
		
		String discretizationStep = "";
		
		for (int i=0; i< ds.length;i++)
			discretizationStep = discretizationStep + String.valueOf(ds[i]) + " ";
		
		properties.setProperty("Variables.discretizationStep", discretizationStep);	  
  }
	
	public static double getDiscretizationStep(int indexOfVariable) {
		int noOfVar = Integer.parseInt(properties.getProperty("Problem.numberOfVariables"));
		double [] steps = new double[noOfVar];
		String stringSteps = properties.getProperty("Variables.discretizationStep");
		String [] stringTableOfSteps = stringSteps.split(" ");
		int counter =0;
		
		for(int i=0; i< stringTableOfSteps.length; i++)
		{
			if(stringTableOfSteps[i] != ""  && stringTableOfSteps[i] != " " && stringTableOfSteps[i].length() > 0)
			{
				steps[counter] = Double.parseDouble(stringTableOfSteps[i]);
				counter++;
			}
		}
		
		return steps[indexOfVariable];
	}

	public static int getNumberOfExactEvaluations() {
		
		return Integer.parseInt(properties.getProperty("DEMONN.exactEvaluation"));
	}

	public static int getNumberOfApproximateEvaluations() {

		return Integer.parseInt(properties.getProperty("DEMONN.approximateEvaluation"));
	}

	public static double[] getLimitsForHypervolume() {

		int noOfObj = Integer.parseInt(properties.getProperty("Problem.numberOfObjectives"));
		double [] limits = new double[noOfObj];
		String stringLimits = properties.getProperty("Problem.hypervolumeCalculationLimits");
		String [] stringTableOfLimits = stringLimits.split(" ");
		int counter =0;
		
		for(int i=0; i< stringTableOfLimits.length; i++)
		{
			if(stringTableOfLimits[i] != ""  && stringTableOfLimits[i] != " " && stringTableOfLimits[i].length() > 0)
			{
				limits[counter] = Double.parseDouble(stringTableOfLimits[i]);
				counter++;
			}
		}
		
		return limits;
		
  }	
	
	public static void setFrontForHypervolume(double[] topBorders) {
	  
		hypervolumeTrueFront =new double[topBorders.length][topBorders.length];
		
		for(int i=0; i< topBorders.length; i++)
			for(int j=0; j< topBorders.length; j++)
				hypervolumeTrueFront[i][j] = 0;

		
		for(int i=0; i< topBorders.length; i++)
			hypervolumeTrueFront[i][i] = topBorders[i];
		
		if(properties.getProperty("Problem.name").contains("Kursawe"))
			hypervolumeTrueFront[0][1] = -31;
  }
	
	public static double [] [] getFrontForHypervolume() {
		
		return hypervolumeTrueFront;
	}
}

