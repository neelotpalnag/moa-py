/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.init;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Time;
import java.util.HashMap;
import java.util.Properties;
import java.util.Scanner;
import jmetal.core.Operator;

import com.mathworks.toolbox.javabuilder.MWException;

import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.metaheuristics.abyss.AbYSS;
import jmetal.metaheuristics.cellde.CellDE;
import jmetal.metaheuristics.demo.DEMO;
import jmetal.metaheuristics.GPDEMO.GPDEMO;
import jmetal.metaheuristics.demoNN.DEMONN;
import jmetal.metaheuristics.demoNN.DEMONNGP;
import jmetal.metaheuristics.demoUncertaintyCompare.DemoUncertaintyCompare;
import jmetal.metaheuristics.exhaustiveSearch.ExhaustiveSearch;
import jmetal.metaheuristics.exhaustiveSearch.ExhaustiveSearchWFG1;
import jmetal.metaheuristics.fastPGA.FastPGA;
import jmetal.metaheuristics.gde3.GDE3;
import jmetal.metaheuristics.ibea.IBEA;
import jmetal.metaheuristics.nsgaII.NSGAII;
import jmetal.metaheuristics.nsgaII.NSGAIIUncertaintyComparison;
import jmetal.metaheuristics.omopso.OMOPSO;
import jmetal.metaheuristics.paes.PAES;
import jmetal.metaheuristics.pesa2.PESA2;
import jmetal.metaheuristics.randomSearch.RandomSearch;
import jmetal.metaheuristics.singleObjective.differentialEvolution.DE;
import jmetal.metaheuristics.singleObjective.evolutionStrategy.ElitistES;
import jmetal.metaheuristics.singleObjective.evolutionStrategy.NonElitistES;
import jmetal.metaheuristics.singleObjective.particleSwarmOptimization.PSO;
import jmetal.metaheuristics.smpso.SMPSO;
import jmetal.metaheuristics.smsemoa.SMSEMOA;
import jmetal.metaheuristics.spea2.SPEA2;
import jmetal.models.ConditionalRF;
import jmetal.models.GP;
import jmetal.models.ModelFactory;
import jmetal.models.NN;
import jmetal.models.RBFNet;
import jmetal.models.RF;
import jmetal.models.SPGPHT;
import jmetal.models.SPGPm;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.localSearch.MutationLocalSearch;
import jmetal.operators.mutation.Mutation;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.mutation.NonUniformMutation;
import jmetal.operators.mutation.UniformMutation;
import jmetal.operators.selection.BinaryTournament;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.BNH;
import jmetal.problems.ConstrEx;
import jmetal.problems.EKGSimulator;
import jmetal.problems.ExternEvaluator;
import jmetal.problems.Fonseca;
import jmetal.problems.Golinski;
import jmetal.problems.IntRealProblem;
import jmetal.problems.Kursawe;
import jmetal.problems.OKA1;
import jmetal.problems.OKA2;
import jmetal.problems.Osyczka;
import jmetal.problems.Osyczka2;
import jmetal.problems.Poloni;
import jmetal.problems.Schaffer;
import jmetal.problems.Srinivas;
import jmetal.problems.StoreSteel;
import jmetal.problems.Tanaka;
import jmetal.problems.Viennet2;
import jmetal.problems.Viennet3;
import jmetal.problems.Viennet4;
import jmetal.problems.Water;
import jmetal.problems.DTLZ.DTLZ1;
import jmetal.problems.DTLZ.DTLZ1a;
import jmetal.problems.DTLZ.DTLZ2;
import jmetal.problems.DTLZ.DTLZ3;
import jmetal.problems.DTLZ.DTLZ4;
import jmetal.problems.DTLZ.DTLZ5;
import jmetal.problems.DTLZ.DTLZ6;
import jmetal.problems.DTLZ.DTLZ7;
import jmetal.problems.LZ09.LZ09_F1;
import jmetal.problems.LZ09.LZ09_F2;
import jmetal.problems.LZ09.LZ09_F3;
import jmetal.problems.LZ09.LZ09_F4;
import jmetal.problems.LZ09.LZ09_F5;
import jmetal.problems.LZ09.LZ09_F6;
import jmetal.problems.LZ09.LZ09_F7;
import jmetal.problems.LZ09.LZ09_F8;
import jmetal.problems.LZ09.LZ09_F9;
import jmetal.problems.WFG.WFG1;
import jmetal.problems.WFG.WFG2;
import jmetal.problems.WFG.WFG3;
import jmetal.problems.WFG.WFG4;
import jmetal.problems.WFG.WFG5;
import jmetal.problems.WFG.WFG6;
import jmetal.problems.WFG.WFG7;
import jmetal.problems.WFG.WFG8;
import jmetal.problems.WFG.WFG9;
import jmetal.problems.ZDT.ZDT1;
import jmetal.problems.ZDT.ZDT2;
import jmetal.problems.ZDT.ZDT3;
import jmetal.problems.ZDT.ZDT4;
import jmetal.problems.ZDT.ZDT5;
import jmetal.problems.ZDT.ZDT6;
import jmetal.problems.cec2009Competition.CEC2009_UF1;
import jmetal.problems.cec2009Competition.CEC2009_UF10;
import jmetal.problems.cec2009Competition.CEC2009_UF2;
import jmetal.problems.cec2009Competition.CEC2009_UF3;
import jmetal.problems.cec2009Competition.CEC2009_UF4;
import jmetal.problems.cec2009Competition.CEC2009_UF5;
import jmetal.problems.cec2009Competition.CEC2009_UF6;
import jmetal.problems.cec2009Competition.CEC2009_UF7;
import jmetal.problems.cec2009Competition.CEC2009_UF8;
import jmetal.problems.cec2009Competition.CEC2009_UF9;
import jmetal.problems.singleObjective.Griewank;
import jmetal.problems.singleObjective.OneMax;
import jmetal.problems.singleObjective.Sphere;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.JMException;
import jmetal.util.comparators.FPGAFitnessComparator;
import jmetal.util.comparators.FitnessComparator;

/**
 *
 * @author Miha Mlakar & Rok Prodan
 */
public class init {
    
    public static Problem createProblem(Properties properties) throws ClassNotFoundException{
        Problem problem=null;
        
        int numberOfVariables=Integer.parseInt(properties.getProperty("Problem.numberOfVariables"));
        String problemName=properties.getProperty("Problem.name");
        
        if (problemName.compareTo("DTLZ1") == 0) {          // 
            problem = new DTLZ1("Real");
        } else if (problemName.compareTo("DTLZ1a") == 0) {          // 
            problem = new DTLZ1a("Real");    
        } else if (problemName.compareTo("DTLZ2") == 0) {    // DTLZ2
            problem = new DTLZ2("Real");
        } else if (problemName.compareTo("DTLZ3") == 0) {    // DTLZ3  
            problem = new DTLZ3("Real");
        } else if (problemName.compareTo("DTLZ4") == 0) {    // DTLZ4
            problem = new DTLZ4("Real");
        } else if (problemName.compareTo("DTLZ5") == 0) {    // DTLZ5
            problem = new DTLZ5("Real");
        } else if (problemName.compareTo("DTLZ6") == 0) {    // DTLZ6
            problem = new DTLZ6("Real");
        } else if (problemName.compareTo("DTLZ7") == 0) {    // DTLZ7
            problem = new DTLZ7("Real", 23, 4);
        } else if (problemName.compareTo("LZ09") == 0) {    // LZ09
            // ni podrazred razreda Problem
            // LZ09 (int nvar, int nobj, int ptype, int dtype, int ltype)
        } else if (problemName.compareTo("LZ09_F1") == 0) {   // LZ09_F1
            problem = new LZ09_F1("Real");
        } else if (problemName.compareTo("LZ09_F2") == 0) {  // LZ09_F2
            problem = new LZ09_F2("Real");
        } else if (problemName.compareTo("LZ09_F3") == 0) {  // LZ09_F3
            problem = new LZ09_F3("Real");
        } else if (problemName.compareTo("LZ09_F4") == 0) {  // LZ09_F4
            problem = new LZ09_F4("Real");
        } else if (problemName.compareTo("LZ09_F5") == 0) {  // LZ09_F5
            problem = new LZ09_F5("Real");
        } else if (problemName.compareTo("LZ09_F6") == 0) {  // LZ09_F6
            problem = new LZ09_F6("Real");
        } else if (problemName.compareTo("LZ09_F7") == 0) {  // LZ09_F7
            problem = new LZ09_F7("Real");
        } else if (problemName.compareTo("LZ09_F8") == 0) {  // LZ09_F8
            problem = new LZ09_F8("Real");
        } else if (problemName.compareTo("LZ09_F9") == 0) {  // LZ09_F9
            problem = new LZ09_F9("Real");
        } else if (problemName.compareTo("Shapes") == 0) {   // Shapes
            // ni podrazred razreda Problem
        } else if (problemName.compareTo("Transformations") == 0) {  // Transformations
            // ni podrazred razreda Problem
        } else if (problemName.compareTo("WFG") == 0) {  // WFG
            // abstract class
        } else if (problemName.compareTo("WFG1") == 0) { // WFG1
            problem = new WFG1("Real");												
        } else if (problemName.compareTo("WFG2") == 0) { // WFG2
            problem = new WFG2("Real");
        } else if (problemName.compareTo("WFG3") == 0) { // WFG3
            problem = new WFG3("Real");
        } else if (problemName.compareTo("WFG4") == 0) { // WFG4
            problem = new WFG4("Real");
        } else if (problemName.compareTo("WFG5") == 0) { // WFG5
            problem = new WFG5("Real");
        } else if (problemName.compareTo("WFG6") == 0) { // WFG6
            problem = new WFG6("Real");
        } else if (problemName.compareTo("WFG7") == 0) { // WFG7
            problem = new WFG7("Real");
        } else if (problemName.compareTo("WFG8") == 0) { // WFG8
            problem = new WFG8("Real");
        } else if (problemName.compareTo("WFG9") == 0) { // WFG9
            problem = new WFG9("Real");
        } else if (problemName.compareTo("ZDT1") == 0) { // ZDT1
            problem = new ZDT1("Real");
        } else if (problemName.compareTo("ZDT2") == 0) { // ZDT2
            problem = new ZDT2("Real");
        } else if (problemName.compareTo("ZDT3") == 0) { // ZDT3
            problem = new ZDT3("Real");
        } else if (problemName.compareTo("ZDT4") == 0) { // ZDT4
            problem = new ZDT4("Real");
        } else if (problemName.compareTo("ZDT5") == 0) { // ZDT5
            problem = new ZDT5("Real");
        } else if (problemName.compareTo("ZDT6") == 0) { // ZDT6
            problem = new ZDT6("Real");
        } else if (problemName.compareTo("CEC2009_UF1") == 0) { // CEC2009_UF1
            problem = new CEC2009_UF1("Real");
        } else if (problemName.compareTo("CEC2009_UF2") == 0) { // CEC2009_UF2
            problem = new CEC2009_UF2("Real");
        } else if (problemName.compareTo("CEC2009_UF3") == 0) { // CEC2009_UF3
            problem = new CEC2009_UF3("Real");
        } else if (problemName.compareTo("CEC2009_UF4") == 0) { // CEC2009_UF4
            problem = new CEC2009_UF4("Real");
        } else if (problemName.compareTo("CEC2009_UF5") == 0) { // CEC2009_UF5
            problem = new CEC2009_UF5("Real");
        } else if (problemName.compareTo("CEC2009_UF6") == 0) { // CEC2009_UF6
            problem = new CEC2009_UF6("Real");
        } else if (problemName.compareTo("CEC2009_UF7") == 0) { // CEC2009_UF7
            problem = new CEC2009_UF7("Real");
        } else if (problemName.compareTo("CEC2009_UF8") == 0) { // CEC2009_UF8
            problem = new CEC2009_UF8("Real");
        } else if (problemName.compareTo("CEC2009_UF9") == 0) { // CEC2009_UF9
            problem = new CEC2009_UF9("Real");
        } else if (problemName.compareTo("CEC2009_UF10") == 0) { // CEC2009_UF10
            problem = new CEC2009_UF10("Real");
        } else if (problemName.compareTo("Griewank") == 0) { // Griewank
            problem = new Griewank("Real", numberOfVariables);
        } else if (problemName.compareTo("OneMax") == 0) { // OneMax
            // OneMax(Integer numberOfBits)
            problem = new OneMax(1);
        } else if (problemName.compareTo("Sphere") == 0) { // Sphere
            problem = new Sphere("Real", numberOfVariables);
        } else if (problemName.compareTo("TSP") == 0) { // TSP
            // TSP(String filename) throws FileNotFoundException, IOException, ClassNotFoundException
            // It accepts data files from TSPLIB (The file containing the definition of the problem)
            // problem = new TSP("");
        } else if (problemName.compareTo("ConstrEx") == 0) { // ConstrEx
            problem = new ConstrEx("Real");
        } else if (problemName.compareTo("Fonseca") == 0) { // Fonseca
            problem = new Fonseca("Real");
        } else if (problemName.compareTo("Golinski") == 0) { // Golinski
            problem = new Golinski("Real");
        } else if (problemName.compareTo("IntRealProblem") == 0) { // IntRealProblem
            problem = new IntRealProblem("IntReal");
        } else if (problemName.compareTo("Kursawe") == 0) { // Kursawe
            problem = new Kursawe("Real", numberOfVariables);
        } else if (problemName.compareTo("OKA1") == 0) { // OKA1
            problem = new OKA1("Real");
        } else if (problemName.compareTo("OKA2") == 0) { // OKA2
            problem = new OKA2("Real");
        } else if (problemName.compareTo("Osyczka") == 0) { // Osyczka2
            problem = new Osyczka("Real");
        } else if (problemName.compareTo("Osyczka2") == 0) { // Osyczka2
          problem = new Osyczka2("Real");
        } else if (problemName.compareTo("Poloni") == 0) { // Poloni
            problem = new Poloni("Real");
        } else if (problemName.compareTo("ProblemFactory") == 0) { // ProblemFactory
            // Class ProblemFactory represents a factory for problems, not a problem.
        } else if (problemName.compareTo("Schaffer") == 0) { // Schaffer
            problem = new Schaffer("Real");
        } else if (problemName.compareTo("Srinivas") == 0) { // Srinivas
            problem = new Srinivas("Real");
        } else if (problemName.compareTo("Tanaka") == 0) { // Tanaka
            problem = new Tanaka("Real");
        } else if (problemName.compareTo("Viennet2") == 0) { // Viennet2
            problem = new Viennet2("Real");
        } else if (problemName.compareTo("Viennet3") == 0) { // Viennet3
            problem = new Viennet3("Real");
        } else if (problemName.compareTo("Viennet4") == 0) { // Viennet4
            problem = new Viennet4("Real");
        } else if (problemName.compareTo("Water") == 0) { // Water
            problem = new Water("Real");
        } else if (problemName.compareTo("BNH") == 0) { // Water
          problem = new BNH("Real");
        } else if (problemName.compareTo("Extern_Evaluator") == 0) { // Extern_Evaluator	
            problem = new ExternEvaluator();
        } else if (problemName.compareTo("StoreSteel") == 0) { // Extern_Evaluator	
          problem = new StoreSteel();
          double [] topBorders = {1, 2, 15};
          Params.setFrontForHypervolume(topBorders);
        } else if (problemName.compareTo("EKGSimulator") == 0) { // Extern_Evaluator	
          problem = new EKGSimulator();
          double [] topBorders = {1, 2, 15};
          Params.setFrontForHypervolume(topBorders);
        } else {
            System.out.printf("Name %s of the PROBLEM in INPUT.txt is invalid", problemName);
        }

        double [] topBorders = new double[problem.getNumberOfObjectives()];
        for(int i=0; i<problem.getNumberOfObjectives(); i++)
        	topBorders[i] = 1000;
        
        if(problem.getName() == "StoreSteel")
        	topBorders = new double[]{1, 2, 15};
        
        else if(problem.getName() == "ExternEvaluator")
        	topBorders = Params.getLimitsForHypervolume();
        
        else if(problem.getName() == "Sphere")
        	topBorders = new double[]{2};
        	
        else if(problem.getName() == "Fonseca")
        	topBorders = new double[]{1, 1};
        
        else if(problem.getName() == "Golinski")
        	topBorders = new double[]{6000, 2000};
        
        else if(problem.getName() == "Kursawe")
        	topBorders = new double[]{-31, 2}; 	//-21, 2
        
        else if(problem.getName() == "Osyczka")
        	topBorders = new double[]{0, 80};
        
        else if(problem.getName() == "Osyczka2")
        	topBorders = new double[]{0, 80};
        
        else if(problem.getName() == "Poloni")
        	topBorders = new double[]{18, 30};
        
        else if(problem.getName() == "Srinivas")
        	topBorders = new double[]{250, 50};
        
        else if(problem.getName() == "Tanaka")
        	topBorders = new double[]{2, 2};
        
        else if(problem.getName() == "Viennet2")
        	topBorders = new double[]{5, -18, -14};
        	
        else if(problem.getName().contains("WFG"))
        {
        	for(int i=0; i<problem.getNumberOfObjectives(); i++)
          	topBorders[i] = 10;
        }
        
        else if(problem.getName().contains("DTLZ3"))
        	topBorders = new double[]{100, 100, 100};
        
        else if(problem.getName().contains("DTLZ"))
        {
        	for(int i=0; i<problem.getNumberOfObjectives(); i++)
          	topBorders[i] = 10;
        }
        
        else if(problem.getName().contains("BNH"))
        	topBorders = new double[]{150, 50};
        
        else if(problem.getName().contains("EKGSimulator"))
        	topBorders = new double[]{2, 2};
        
        
        Params.setFrontForHypervolume(topBorders);
        
        //if (problemName.compareTo("Extern_Evaluator") != 0)
        //	problem.setNumberOfFeatures(problem.getNumberOfObjectives());
        
        return problem;
    }
    
    public static Algorithm createAlgorithm(Properties properties, Problem problem) throws Exception {
    	Operator crossover = null;         // Crossover operator
        Operator mutation = null;         // Mutation operator
        Operator selection = null;         // Selection operator

        Operator parentsSelection = null;       // ParentSelection operator for MOCHC algorithm
        Operator newGenerationSelection = null; // NewGenerationSelection for MOCHC algorithm

        Operator improvement = null; // Improvement operator for Abyss algorithm

        Mutation uniformMutation = null;   // UniformMutation for OMOPSO algorithm
        Mutation nonUniformMutation = null;// NonUniFormMutation for OMOPSO algorithm

        Algorithm algorithm = null;
        HashMap parameters;

        String algorithmName = properties.getProperty("Algorithm.name");

        if (algorithmName.compareTo("AbYSS") == 0) {
            // STEP 4. Specify and configure the crossover operator, used in the
            //         solution combination method of the scatter search
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            // STEP 5. Specify and configure the improvement method. We use by default
            //         a polynomial mutation in this method.
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            parameters = new HashMap();
            parameters.put("improvementRounds", Integer.valueOf(properties.getProperty("AbYSS.improvementRounds")));
            parameters.put("problem", problem);
            parameters.put("mutation", mutation);
            improvement = new MutationLocalSearch(parameters); // Operator for improvement

            algorithm = new AbYSS(problem);
        } else if (algorithmName.compareTo("CellDE") == 0) {
            // Crossover operator 
            parameters = new HashMap();
            parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

            parameters = null;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament", parameters);

            algorithm = new CellDE(problem);
        } else if (algorithmName.compareTo("DEMO") == 0) {
            // Crossover operator 
            parameters = new HashMap();
            parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

            parameters = null;
            selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters);

            algorithm = new DEMO(problem, properties);

        } else if (algorithmName.compareTo("DemoUncertaintyCompare") == 0) {
          // Crossover operator 
        	parameters = new HashMap();
		      parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
		      parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
		      crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
		
		      parameters = null;
		      selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters);

		      //set operators_
		      //algorithm.addOperator("crossover",crossover);
	        //algorithm.addOperator("mutation",mutation);
	        //algorithm.addOperator("selection",selection);		 
		      
		      Params.setSelectionOperator(selection);
		      Params.setCrossoverOperator(crossover);
		      
	      	//algorithm = new GPDEMO(problem, ModelFactory.getModel(properties.getProperty("Model.name"), problem.getNumberOfVariables(), problem.getNumberOfObjectives()), properties);
	      	
	      	Model model;
	      	Model modelLocal;
	      	int activeSet;
	      	int activeSetLocal;
	      	int learningType = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	int learningTypeLocal = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	boolean reset = false;
	      	boolean resetLocal = false;
	      	String covfunc;
	      	String covfuncLocal;
	      	double[][] hyp;
	      	double[][] hypLocal;
	      	String val;
	      	int maxOptIter;
	      	double optThreshold;
	      	int maxOptIterLocal;
	      	double optThresholdLocal;
	      	
	      	
	      	// global model
	      	try {
	      		activeSet = Integer.parseInt(properties.getProperty("Model.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSet = -1;
    			}
	      	
	      	try {
	      		learningType = Integer.parseInt(properties.getProperty("Model.learningType"));
	      		reset = (learningType == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("Model.learningType in not set correctly.");
    				System.exit(0);
    			}
	      	
	      	try {
	      		covfunc = properties.getProperty("Model.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfunc = "covSEard";
    			}

	      	try {
	      		val = properties.getProperty("Model.hyp");
	      		if (val != null)
	      			hyp = parseHypVal(val);
	      		else
	      			hyp = null;
	      	}
    			catch (NumberFormatException e) {
    				hyp = null;
    			}

	      	try {
	      		maxOptIter = Integer.parseInt(properties.getProperty("Model.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIter = -500;
    			}

	      	try {
	      		optThreshold = Double.parseDouble(properties.getProperty("Model.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThreshold = 0.0000001;
    			}
	      	
	      	switch (properties.getProperty("Model.name")) {
	      		case "RF":
	      			model = new RF(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), Integer.parseInt(properties.getProperty("RF.onOfTrees")));
	      			break;
	      		case "ConditionalRF":
	      			model = new ConditionalRF(problem.getNumberOfVariables(), problem.getNumberOfObjectives());
	      			break;
	      		case "RBF":
	      			model = new RBFNet(problem.getNumberOfVariables(), problem.getNumberOfObjectives());
	      			break;
	      		case "GP":
	      			model = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet, reset, covfunc, hyp);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPm)model).setResetWin(true);
	      			break;
	      		case "SPGPHT":
	      			if (activeSet > 0)
	      				model = new SPGPHT(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			((SPGPHT)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPHT)model).setHyperparameters(hyp);
	      			((SPGPHT)model).setMaxOptIter(maxOptIter);
	      			//((SPGPHT)model).setOptThreshold(optThreshold);
	      			((SPGPHT)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPHT)model).setResetWin(true);
	      			break;
	      		default:
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), 300);
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPm)model).setResetWin(true);
	      			break;
	      	}

	      	// local model
	      	try {
	      		activeSetLocal = Integer.parseInt(properties.getProperty("ModelLocal.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSetLocal = -1;
    			}
	      	try {
	      		learningTypeLocal = Integer.parseInt(properties.getProperty("ModelLocal.learningType"));
	      		resetLocal = (learningTypeLocal == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("ModelLocal.learningType in not set correctly.");
    				System.exit(0);
    			}	      	
	      	try {
	      		covfuncLocal = properties.getProperty("ModelLocal.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfuncLocal = "covSEard";
    			}
	      	try {
	      		val = properties.getProperty("ModelLocal.hyp");
	      		if (val != null)
	      			hypLocal = parseHypVal(val);
	      		else
	      			hypLocal = null;
	      	}
    			catch (NumberFormatException e) {
    				hypLocal = null;
    			}
	      	try {
	      		maxOptIterLocal = Integer.parseInt(properties.getProperty("ModelLocal.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIterLocal = -500;
    			}
	      	try {
	      		optThresholdLocal = Double.parseDouble(properties.getProperty("ModelLocal.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThresholdLocal = 0.0000001;
    			}
	      	switch (properties.getProperty("ModelLocal.name")) {
	      		case "GP":
	      			modelLocal = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal, resetLocal, covfuncLocal, hypLocal);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0) {
	      				modelLocal = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal);
		      			((SPGPm)modelLocal).setResetHyp(resetLocal);
		      			if (hypLocal != null)
		      				((SPGPm)modelLocal).setHyperparameters(hypLocal);
		      			((SPGPm)modelLocal).setMaxOptIter(maxOptIterLocal);
		      			//((SPGPm)modelLocal).setOptThreshold(optThresholdLocal);
	      			}	      			      		      		
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			break;
	      			
	      		case "RF":
	      			modelLocal =new RF(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), Integer.parseInt(properties.getProperty("RF.onOfTrees")));	    
	      			break;
	      			
	      		default:
	      			modelLocal = null;
	      	}
		      
	      	algorithm = new DemoUncertaintyCompare(problem, model, modelLocal, properties);	          
          
        } else if (algorithmName.compareTo("GPDEMO") == 0) {
		      // Crossover operator 
		      parameters = new HashMap();
		      parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
		      parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
		      crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
		
		      parameters = null;
		      selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters);

		      //set operators_
		      //algorithm.addOperator("crossover",crossover);
	        //algorithm.addOperator("mutation",mutation);
	        //algorithm.addOperator("selection",selection);		 
		      
		      Params.setSelectionOperator(selection);
		      Params.setCrossoverOperator(crossover);
		      
	      	//algorithm = new GPDEMO(problem, ModelFactory.getModel(properties.getProperty("Model.name"), problem.getNumberOfVariables(), problem.getNumberOfObjectives()), properties);
	      	
	      	Model model;
	      	Model modelLocal;
	      	int activeSet;
	      	int activeSetLocal;
	      	int learningType = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	int learningTypeLocal = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	boolean reset = false;
	      	boolean resetLocal = false;
	      	String covfunc;
	      	String covfuncLocal;
	      	double[][] hyp;
	      	double[][] hypLocal;
	      	String val;
	      	int maxOptIter;
	      	double optThreshold;
	      	int maxOptIterLocal;
	      	double optThresholdLocal;
	      	
	      	
	      	// global model
	      	try {
	      		activeSet = Integer.parseInt(properties.getProperty("Model.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSet = -1;
    			}
	      	
	      	try {
	      		learningType = Integer.parseInt(properties.getProperty("Model.learningType"));
	      		reset = (learningType == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("Model.learningType in not set correctly.");
    				System.exit(0);
    			}
	      	
	      	try {
	      		covfunc = properties.getProperty("Model.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfunc = "covSEard";
    			}

	      	try {
	      		val = properties.getProperty("Model.hyp");
	      		if (val != null)
	      			hyp = parseHypVal(val);
	      		else
	      			hyp = null;
	      	}
    			catch (NumberFormatException e) {
    				hyp = null;
    			}

	      	try {
	      		maxOptIter = Integer.parseInt(properties.getProperty("Model.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIter = -500;
    			}

	      	try {
	      		optThreshold = Double.parseDouble(properties.getProperty("Model.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThreshold = 0.0000001;
    			}
	      	
	      	switch (properties.getProperty("Model.name")) {
	      		case "RF":
	      			model = new RF(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), Integer.parseInt(properties.getProperty("RF.onOfTrees")));
	      			break;
	      		case "ConditionalRF":
	      			model = new ConditionalRF(problem.getNumberOfVariables(), problem.getNumberOfObjectives());
	      			break;
	      		case "RBF":
	      			model = new RBFNet(problem.getNumberOfVariables(), problem.getNumberOfObjectives());
	      			break;
	      		case "GP":
	      			model = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet, reset, covfunc, hyp);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPm)model).setResetWin(true);
	      			break;
	      		case "SPGPHT":
	      			if (activeSet > 0)
	      				model = new SPGPHT(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			((SPGPHT)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPHT)model).setHyperparameters(hyp);
	      			((SPGPHT)model).setMaxOptIter(maxOptIter);
	      			//((SPGPHT)model).setOptThreshold(optThreshold);
	      			((SPGPHT)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPHT)model).setResetWin(true);
	      			break;
	      		default:
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), 300);
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
	      			//((SPGPm)model).setResetWin(true);
	      			break;
	      	}

	      	// local model
	      	try {
	      		activeSetLocal = Integer.parseInt(properties.getProperty("ModelLocal.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSetLocal = -1;
    			}
	      	try {
	      		learningTypeLocal = Integer.parseInt(properties.getProperty("ModelLocal.learningType"));
	      		resetLocal = (learningTypeLocal == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("ModelLocal.learningType in not set correctly.");
    				System.exit(0);
    			}	      	
	      	try {
	      		covfuncLocal = properties.getProperty("ModelLocal.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfuncLocal = "covSEard";
    			}
	      	try {
	      		val = properties.getProperty("ModelLocal.hyp");
	      		if (val != null)
	      			hypLocal = parseHypVal(val);
	      		else
	      			hypLocal = null;
	      	}
    			catch (NumberFormatException e) {
    				hypLocal = null;
    			}
	      	try {
	      		maxOptIterLocal = Integer.parseInt(properties.getProperty("ModelLocal.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIterLocal = -500;
    			}
	      	try {
	      		optThresholdLocal = Double.parseDouble(properties.getProperty("ModelLocal.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThresholdLocal = 0.0000001;
    			}
	      	switch (properties.getProperty("ModelLocal.name")) {
	      		case "GP":
	      			modelLocal = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal, resetLocal, covfuncLocal, hypLocal);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0) {
	      				modelLocal = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal);
		      			((SPGPm)modelLocal).setResetHyp(resetLocal);
		      			if (hypLocal != null)
		      				((SPGPm)modelLocal).setHyperparameters(hypLocal);
		      			((SPGPm)modelLocal).setMaxOptIter(maxOptIterLocal);
		      			//((SPGPm)modelLocal).setOptThreshold(optThresholdLocal);
	      			}	      			      		      		
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			break;
	      			
	      		case "RF":
	      			modelLocal =new RF(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), Integer.parseInt(properties.getProperty("RF.onOfTrees")));	    
	      			break;
	      			
	      		default:
	      			modelLocal = null;
	      	}
		      
	      	algorithm = new GPDEMO(problem, model, modelLocal, properties);		      	
		    
        } else if (algorithmName.compareTo("DEMONNGP") == 0) {
		      // Crossover operator 
		      parameters = new HashMap();
		      parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
		      parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
		      crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);
		
		      parameters = null;
		      selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters);

		      //set operators_
		      //algorithm.addOperator("crossover",crossover);
	        //algorithm.addOperator("mutation",mutation);
	        //algorithm.addOperator("selection",selection);		 
		      
		      Params.setSelectionOperator(selection);
		      Params.setCrossoverOperator(crossover);
		      
	      	//algorithm = new GPDEMO(problem, ModelFactory.getModel(properties.getProperty("Model.name"), problem.getNumberOfVariables(), problem.getNumberOfObjectives()), properties);
	      	
	      	Model model;
	      	Model modelLocal;
	      	int activeSet;
	      	int activeSetLocal;
	      	int learningType = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	int learningTypeLocal = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
	      	boolean reset = false;
	      	boolean resetLocal = false;
	      	String covfunc;
	      	String covfuncLocal;
	      	double[][] hyp;
	      	double[][] hypLocal;
	      	String val;
	      	int maxOptIter;
	      	double optThreshold;
	      	int maxOptIterLocal;
	      	double optThresholdLocal;
	      	
	      	
	      	// global model
	      	try {
	      		activeSet = Integer.parseInt(properties.getProperty("Model.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSet = -1;
    			}
	      	
	      	try {
	      		learningType = Integer.parseInt(properties.getProperty("Model.learningType"));
	      		reset = (learningType == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("Model.learningType in not set correctly.");
    				System.exit(0);
    			}
	      	
	      	try {
	      		covfunc = properties.getProperty("Model.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfunc = "covSEard";
    			}

	      	try {
	      		val = properties.getProperty("Model.hyp");
	      		if (val != null)
	      			hyp = parseHypVal(val);
	      		else
	      			hyp = null;
	      	}
    			catch (NumberFormatException e) {
    				hyp = null;
    			}

	      	try {
	      		maxOptIter = Integer.parseInt(properties.getProperty("Model.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIter = -500;
    			}

	      	try {
	      		optThreshold = Double.parseDouble(properties.getProperty("Model.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThreshold = 0.0000001;
    			}
	      	
	      	switch (properties.getProperty("Model.name")) {
	      		case "GP":
	      			model = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet, reset, covfunc, hyp);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			break;
	      		default:
	      			if (activeSet > 0)
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
	      			else
	      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), 300);
	      			((SPGPm)model).setResetHyp(reset);
	      			if (hyp != null)
	      				((SPGPm)model).setHyperparameters(hyp);
	      			((SPGPm)model).setMaxOptIter(maxOptIter);
	      			//((SPGPm)model).setOptThreshold(optThreshold);
	      			break;
	      	}

	      	// local model
	      	try {
	      		activeSetLocal = Integer.parseInt(properties.getProperty("ModelLocal.activeSet"));
	      	}
    			catch (NumberFormatException e) {
    				activeSetLocal = -1;
    			}
	      	try {
	      		learningTypeLocal = Integer.parseInt(properties.getProperty("ModelLocal.learningType"));
	      		resetLocal = (learningTypeLocal == 0) ? true : false;
	      	}
	      	catch (Exception e) {
	      		System.out.println("ModelLocal.learningType in not set correctly.");
    				System.exit(0);
    			}	      	
	      	try {
	      		covfuncLocal = properties.getProperty("ModelLocal.covfunc");
	      	}
    			catch (NumberFormatException e) {
    				covfuncLocal = "covSEard";
    			}
	      	try {
	      		val = properties.getProperty("ModelLocal.hyp");
	      		if (val != null)
	      			hypLocal = parseHypVal(val);
	      		else
	      			hypLocal = null;
	      	}
    			catch (NumberFormatException e) {
    				hypLocal = null;
    			}
	      	try {
	      		maxOptIterLocal = Integer.parseInt(properties.getProperty("ModelLocal.maxOptIter"));
	      	}
    			catch (NumberFormatException e) {
    				maxOptIterLocal = -500;
    			}
	      	try {
	      		optThresholdLocal = Double.parseDouble(properties.getProperty("ModelLocal.optThreshold"));
	      	}
    			catch (NumberFormatException e) {
    				optThresholdLocal = 0.0000001;
    			}
	      	switch (properties.getProperty("ModelLocal.name")) {
	      		case "GP":
	      			modelLocal = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal, resetLocal, covfuncLocal, hypLocal);
	      			break;
	      		case "SPGP":
	      			if (activeSet > 0) {
	      				modelLocal = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSetLocal);
		      			((SPGPm)modelLocal).setResetHyp(resetLocal);
		      			if (hypLocal != null)
		      				((SPGPm)modelLocal).setHyperparameters(hypLocal);
		      			((SPGPm)modelLocal).setMaxOptIter(maxOptIterLocal);
		      			//((SPGPm)modelLocal).setOptThreshold(optThresholdLocal);
	      			}
	      			else
	      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
	      			break;
	      		default:
	      			modelLocal = null;
	      	}
		      
	      	algorithm = new DEMONNGP(problem, model, modelLocal, properties);		      	
		    
	      	
		    } else if(algorithmName.compareTo("DEMONN")==0){
        	parameters = new HashMap() ;
            parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability"))) ;
            parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);   

            // Add the operators to the algorithm
            parameters = null ;
            selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters) ;
        	
            algorithm=new DEMONN(problem, properties);
            
        } else if(algorithmName.compareTo("DENSEA")==0){  
 
        } else if(algorithmName.compareTo("ExhaustiveSearch")==0){  
            algorithm=new ExhaustiveSearch(problem, properties);
        } else if(algorithmName.compareTo("ExhaustiveSearchWFG1")==0){  
          algorithm=new ExhaustiveSearchWFG1(problem, properties);
        } else if(algorithmName.compareTo("FastPGA")==0){  
        	// Mutation and Crossover for Real codification 
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);
            // Mutation and Crossover for Binary codification

            parameters = new HashMap();
            parameters.put("comparator", new FPGAFitnessComparator());
            selection = new BinaryTournament(parameters);

            algorithm = new FastPGA(problem);
        } else if(algorithmName.compareTo("GDE3")==0){ 
            // Crossover operator 
            parameters = new HashMap() ;
            parameters.put("CR", 0.5) ;
            parameters.put("F", 0.5) ;
            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);                   

            // Add the operators to the algorithm
            parameters = null ;
            selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters) ;
            
            algorithm=new GDE3(problem, properties);
            
        } else if(algorithmName.compareTo("IBEA")==0){
            // Mutation and Crossover for Real codification 
            parameters = new HashMap() ;
            parameters.put("probability", 0.9) ;
            parameters.put("distributionIndex", 20.0) ;
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);                   

            parameters = new HashMap() ;
            parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
            parameters.put("distributionIndex", 20.0) ;
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);         

            /* Selection Operator */
            parameters = new HashMap() ; 
            parameters.put("comparator", new FitnessComparator()) ;
            selection = new BinaryTournament(parameters);
            
            algorithm=new IBEA(problem, properties);
            
        } else if(algorithmName.compareTo("MOCHC")==0){  
        	 // Crossover operator
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            crossover = CrossoverFactory.getCrossoverOperator("HUXCrossover", parameters);

            parameters = null;
            parentsSelection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);

            parameters = new HashMap();
            parameters.put("problem", problem);
            newGenerationSelection = SelectionFactory.getSelectionOperator("RankingAndCrowdingSelection", parameters);

            // Mutation operator
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")));
            mutation = MutationFactory.getMutationOperator("BitFlipMutation", parameters);
        } else if(algorithmName.compareTo("NSGAII")==0){ 
             // Mutation and Crossover for Real codification 
            parameters = new HashMap() ;
            parameters.put("probability", 0.9) ;
            parameters.put("distributionIndex", 20.0) ;
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);                   

            parameters = new HashMap() ;
            parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
            parameters.put("distributionIndex", 20.0) ;
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    

            // Selection Operator 
            parameters = null ;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters); 
            
            algorithm=new NSGAII(problem, properties);
            
        } else if(algorithmName.compareTo("NSGAIIUncertaintyComparison")==0){ 
		          // Mutation and Crossover for Real codification 
		         parameters = new HashMap() ;
		         parameters.put("probability", 0.9) ;
		         parameters.put("distributionIndex", 20.0) ;
		         crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);                   
		
		         parameters = new HashMap() ;
		         parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
		         parameters.put("distributionIndex", 20.0) ;
		         mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    
		
		         // Selection Operator 
		         parameters = null ;
		         selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters); 
		         
		         Model model;
			      	Model modelLocal;
			      	int activeSet;
			      	int activeSetLocal;
			      	int learningType = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
			      	int learningTypeLocal = 0; //0 = learn from the begining  1 = update previously learned hiperparameters
			      	boolean reset = false;
			      	boolean resetLocal = false;
			      	String covfunc;
			      	String covfuncLocal;
			      	double[][] hyp;
			      	double[][] hypLocal;
			      	String val;
			      	int maxOptIter;
			      	double optThreshold;
			      	int maxOptIterLocal;
			      	double optThresholdLocal;
			      	
			      	
			      	// global model
			      	try {
			      		activeSet = Integer.parseInt(properties.getProperty("Model.activeSet"));
			      	}
		    			catch (NumberFormatException e) {
		    				activeSet = -1;
		    			}
			      	
			      	try {
			      		learningType = Integer.parseInt(properties.getProperty("Model.learningType"));
			      		reset = (learningType == 0) ? true : false;
			      	}
			      	catch (Exception e) {
			      		System.out.println("Model.learningType in not set correctly.");
		    				System.exit(0);
		    			}
			      	
			      	try {
			      		covfunc = properties.getProperty("Model.covfunc");
			      	}
		    			catch (NumberFormatException e) {
		    				covfunc = "covSEard";
		    			}

			      	try {
			      		val = properties.getProperty("Model.hyp");
			      		if (val != null)
			      			hyp = parseHypVal(val);
			      		else
			      			hyp = null;
			      	}
		    			catch (NumberFormatException e) {
		    				hyp = null;
		    			}

			      	try {
			      		maxOptIter = Integer.parseInt(properties.getProperty("Model.maxOptIter"));
			      	}
		    			catch (NumberFormatException e) {
		    				maxOptIter = -500;
		    			}

			      	try {
			      		optThreshold = Double.parseDouble(properties.getProperty("Model.optThreshold"));
			      	}
		    			catch (NumberFormatException e) {
		    				optThreshold = 0.0000001;
		    			}
			      	
		         switch (properties.getProperty("Model.name")) {
			      		case "RF":
			      			model = new RF(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), Integer.parseInt(properties.getProperty("RF.onOfTrees")));
			      			break;			      		
			      		case "GP":
			      			model = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet, reset, covfunc, hyp);
			      			break;
			      		case "SPGP":
			      			if (activeSet > 0)
			      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
			      			else
			      				throw new Exception("Property 'ActiveSet' limit is obligatory for a SPGP model.");
			      			((SPGPm)model).setResetHyp(reset);
			      			if (hyp != null)
			      				((SPGPm)model).setHyperparameters(hyp);
			      			((SPGPm)model).setMaxOptIter(maxOptIter);
			      			//((SPGPm)model).setOptThreshold(optThreshold);
			      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
			      			//((SPGPm)model).setResetWin(true);
			      			break;
			      		default:
			      			if (activeSet > 0)
			      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), activeSet);
			      			else
			      				model = new SPGPm(problem.getNumberOfVariables(), problem.getNumberOfObjectives(), 300);
			      			((SPGPm)model).setResetHyp(reset);
			      			if (hyp != null)
			      				((SPGPm)model).setHyperparameters(hyp);
			      			((SPGPm)model).setMaxOptIter(maxOptIter);
			      			//((SPGPm)model).setOptThreshold(optThreshold);
			      			((SPGPm)model).setWindowSize(Integer.valueOf(properties.getProperty("Model.windowSize")));
			      			//((SPGPm)model).setResetWin(true);
			      			break;
			      	}
		         
		         model.setName(properties.getProperty("Model.name"));
		         algorithm=new NSGAIIUncertaintyComparison(problem, properties, model);
		         
		   }   else if (algorithmName.compareTo("OMOPSO") == 0) {
            Integer maxIterations = Integer.valueOf(properties.getProperty("OMOPSO.maxIterations"));
            Double perturbationIndex = Double.valueOf(properties.getProperty("OMOPSO.perturbationIndex"));
            Double mutationProbability = Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables();

            parameters = new HashMap();
            parameters.put("probability", mutationProbability);
            parameters.put("perturbation", perturbationIndex);
            uniformMutation = new UniformMutation(parameters);

            parameters = new HashMap();
            parameters.put("probability", mutationProbability);
            parameters.put("perturbation", perturbationIndex);
            parameters.put("maxIterations", maxIterations);
            nonUniformMutation = new NonUniformMutation(parameters);

            algorithm = new OMOPSO(problem);
        } else if (algorithmName.compareTo("PAES") == 0) {
            // Mutation (Real variables)
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            algorithm = new PAES(problem);
        } else if (algorithmName.compareTo("PESA2") == 0) {
            // Mutation and Crossover for Real codification 
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            // Mutation and Crossover Binary codification
            /*
            crossover = CrossoverFactory.getCrossoverOperator("SinglePointCrossover");                   
            crossover.setParameter("probability",0.9);                   
            mutation = MutationFactory.getMutationOperator("BitFlipMutation");                    
            mutation.setParameter("probability",1.0/80);
             */

            algorithm = new PESA2(problem);
        } else if (algorithmName.compareTo("RandomSearch") == 0) {
            algorithm = new RandomSearch(problem, properties);
        } else if (algorithmName.compareTo("DE") == 0) {
            // Crossover operator 
            parameters = new HashMap();
            parameters.put("CR", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("F", Double.valueOf(properties.getProperty("DEMO.weight")));
            parameters.put("DE_VARIANT", "rand/1/bin");

            crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

            parameters = null;
            selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters);

            algorithm = new DE(problem);   // Asynchronous cGA
        } else if (algorithmName.compareTo("ElitistES") == 0) {
            int bits; // Length of bit string in the OneMax problem
            bits = Integer.valueOf(properties.getProperty("ElitismES.bits"));

            int mu;
            int lambda;
            // Requirement: lambda must be divisible by mu
            mu = Integer.valueOf(properties.getProperty("ElitismES.mu"));
            lambda = Integer.valueOf(properties.getProperty("ElitismES.lambda"));

            /* Mutation and Crossover for Real codification */
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / bits);
            mutation = MutationFactory.getMutationOperator("BitFlipMutation", parameters);

            algorithm = new ElitistES(problem, mu, lambda);
        } else if (algorithmName.compareTo("NonElitistES") == 0) {
            int bits; // Length of bit string in the OneMax problem
            bits = Integer.valueOf(properties.getProperty("ElitismES.bits"));

            int mu;
            int lambda;
            // Requirement: lambda must be divisible by mu
            mu = Integer.valueOf(properties.getProperty("ElitismES.mu"));
            lambda = Integer.valueOf(properties.getProperty("ElitismES.lambda"));

            /* Mutation and Crossover for Real codification */
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / bits);
            mutation = MutationFactory.getMutationOperator("BitFlipMutation", parameters);

            algorithm = new NonElitistES(problem, mu, lambda);
        } else if (algorithmName.compareTo("PSO") == 0) {
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            algorithm = new PSO(problem);
        } else if (algorithmName.compareTo("SMPSO") == 0) {
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            algorithm = new SMPSO(problem);
        } else if (algorithmName.compareTo("SMSEMOA") == 0) {
            // Mutation and Crossover for Real codification 
            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.crossoverProbability")));
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

            parameters = new HashMap();
            parameters.put("probability", Double.valueOf(properties.getProperty("Algorithm.mutationProbability")) / problem.getNumberOfVariables());
            parameters.put("distributionIndex", Double.valueOf(properties.getProperty("Algorithm.distributionIndex")));
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

            // Selection Operator
            parameters = null;
            selection = SelectionFactory.getSelectionOperator("RandomSelection", parameters);
            // also possible
            // selection = SelectionFactory.getSelectionOperator("BinaryTournament2");

            algorithm = new SMSEMOA(problem);
        } else if(algorithmName.compareTo("SPEA2")==0){ 
             // Mutation and Crossover for Real codification 
            parameters = new HashMap() ;
            parameters.put("probability", 0.9) ;
            parameters.put("distributionIndex", 20.0) ;
            crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);                   

            parameters = new HashMap() ;
            parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
            parameters.put("distributionIndex", 20.0) ;
            mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    

            // Selection operator 
            parameters = null ;
            selection = SelectionFactory.getSelectionOperator("BinaryTournament", parameters);
            
            algorithm=new SPEA2(problem, properties);
            
        } else{
            System.out.printf("Name %s of the ALGORITHM in INPUT.txt is invalid", algorithmName);
        }
        
        // Add the operators to the algorithm
        algorithm.addOperator("crossover",crossover);
        algorithm.addOperator("mutation",mutation);
        algorithm.addOperator("selection",selection);
        
        if (parentsSelection != null && newGenerationSelection != null) {                // MOCHC algorithm
            algorithm.addOperator("parentsSelection", parentsSelection);
            algorithm.addOperator("newGenerationSelection", newGenerationSelection);
        } else if (improvement != null) {                                                // Abyss algorithm
            algorithm.addOperator("improvement", improvement);
        } else if (uniformMutation != null && nonUniformMutation != null) {              // OMOPSO algorithm
            algorithm.addOperator("uniformMutation", uniformMutation);
            algorithm.addOperator("nonUniformMutation", nonUniformMutation);
        }
        
        return algorithm;
    }
            
    private static double[][] parseHypVal(String hyps) throws Exception {
    	String[] hypsarr;
    	String[] hypsval;
    	double[][] hyp;
    	int len;
    	
    	len = -1;
    	hypsarr = hyps.split(" ");
    	hyp = new double[0][0];
    	for (int i = 0; i < hypsarr.length; i++) {
    		hypsval = hypsarr[i].split(";");
    		if (len == -1) {
    			len = hypsval.length;
    			hyp = new double[len][hypsarr.length];
    		}
    		if (len != hypsval.length)
    			throw new Exception("Hiperparametri se ne ujemajo!");
    		for (int j = 0; j < hypsval.length; j++) {
    			hyp[j][i] = Integer.parseInt(hypsval[j]);
    		}
    	}
	    return hyp;
    }

		public static void main(String[] args) throws Exception{
            	
        Problem   problem   ;         // The problem to solve
        Algorithm algorithm ;         // The algorithm to use
        String propertyName;
        String propertyValue = "";
               
        Properties properties;
        properties =new Properties();
        String line=null;
        Scanner s=null;
        BufferedReader br=null;
        SolutionSet population=new SolutionSet();
        
        QualityIndicator indicators=null; // Object to get quality indicators
        
        try{
            br=new BufferedReader(new FileReader(args[0]));    
            while((line=br.readLine())!=null){
                if(line.equals("") || line.charAt(0)=='@'){
                    //Comments are ignored
                }else{
                    s=new Scanner(line);
                    while(s.hasNext()){
                    	propertyName = s.next();
                    	if(propertyName.equals("Variables.lowerLimit"))
                    	{
                    		propertyValue = "";
                    		while(s.hasNext())
                    			propertyValue =propertyValue + " " + s.next();
                    		properties.setProperty(propertyName, propertyValue);
                    	}
                    	else if (propertyName.equals("Variables.upperLimit"))
                    	{
                    		propertyValue = "";
                    		while(s.hasNext())
                    			propertyValue =propertyValue + " " + s.next();
                    		properties.setProperty(propertyName, propertyValue);
                    	}
                    	else if (propertyName.equals("Variables.discretizationStep"))
                    	{
                    		propertyValue = "";
                    		while(s.hasNext())
                    			propertyValue =propertyValue + " " + s.next();
                    		properties.setProperty(propertyName, propertyValue);
                    	}
                    	else if (propertyName.equals("Problem.hypervolumeCalculationLimits"))
                    	{
                    		propertyValue = "";
                    		while(s.hasNext())
                    			propertyValue =propertyValue + " " + s.next();
                    		properties.setProperty(propertyName, propertyValue);
                    	}
                    	// GP models hyperparameter values
                    	else if(propertyName.equals("Model.hyp") || propertyName.equals("ModelLocal.hyp"))
                    	{
                    		propertyValue = "";
                    		while(s.hasNext())
                    			propertyValue += s.next() + " ";
                    		propertyValue.substring(0,  propertyValue.length()-1);
                    		properties.setProperty(propertyName, propertyValue);
                    	}
                    	else
                    		properties.setProperty(propertyName, s.next());
                    }
                }
            }
        }
        finally{
            if(br!=null){
                br.close();               
            }
            if(s!=null){
                s.close();
            }
        }
        
        if(args.length > 1){
            if ((args.length % 2) == 1){
                for(int i=1; i<args.length; i+=2){
                    properties.setProperty(args[i], args[i+1]); //Overriding properties
                }
            }else{
                System.out.println("INCORRECT USE!");
                System.out.println("java init <INPUT.txt> <parameter_1> <value_1> ... <parameter_N> <value_N>");
                System.out.println("See INPUT.txt for names of parameters");
            }
        }
        
        Params.setProperties(properties);
        
        problem=createProblem(properties);
        algorithm=createAlgorithm(properties,problem);
        
        //repairing the problem and algorithm parameters for the build in problems
        properties.setProperty("Problem.numberOfVariables", String.valueOf(problem.getNumberOfVariables()));
        properties.setProperty("Problem.numberOfObjectives", String.valueOf(problem.getNumberOfObjectives()));
        properties.setProperty("Problem.numberOfFeatures", String.valueOf(problem.getNumberOfFeatures()));               
        
        // Add the indicator object to the algorithm
        algorithm.setInputParameter("indicators", indicators);
                
        /*
         * Write parameters to Log, Front and Generation file
         */
        String logFileName = properties.getProperty("Output.logFileName");
        int logMode=Integer.parseInt(properties.getProperty("Output.logMode"));
        
        if(logMode==1){ //Overwrite - log file
            population.printParametersToFile(logFileName, properties, false);
        }else if(logMode==2){ //Append - log file
            population.printParametersToFile(logFileName, properties, true);
        }
        
        String frontFileName = properties.getProperty("Output.frontFileName");
        int frontMode=Integer.parseInt(properties.getProperty("Output.frontMode"));
        
        if(frontMode==1){ //Overwrite 
            population.printParametersToFile(frontFileName, properties, false);
        }else if(frontMode==2){ //Append 
            population.printParametersToFile(frontFileName, properties, true);
        }
        
        String genFileName = properties.getProperty("Output.genFileName");
        int genMode=Integer.parseInt(properties.getProperty("Output.genMode"));
        
        if(genMode==1){ //Overwrite 
            population.printParametersToFile(genFileName, properties, false);
            population.printToFile(genFileName, "#gen_number	num_feasible	num_unfeasible	num_non_evaluated	avg_violation	num_nondominated	exact_hyp	spread	c_better	p_better	incomparable	additionalEvaluations	additionalEvaluationsDuringEvolution");
        }else if(genMode==2){ //Append 
            population.printParametersToFile(genFileName, properties, true);
            population.printToFile(genFileName, "#gen_number	num_feasible	num_unfeasible	num_non_evaluated	avg_violation	num_nondominated	exact_hyp	spread	c_better	p_better	incomparable 	additionalEvaluations	additionalEvaluationsDuringEvolution");
        }
        
        
        /*
         * Execute the Algorithm
         */
        long initTime = System.currentTimeMillis();
        Params.setStartTime(initTime);
        System.out.println("Zaèetek = " + new Time(initTime));
        population = algorithm.execute(); 
        long endTime = System.currentTimeMillis();
        System.out.println("Konec = " + new Time(endTime));
        long diff = endTime - initTime;
        
        long secondInMillis = 1000;
        long minuteInMillis = secondInMillis * 60;
        long hourInMillis = minuteInMillis * 60;
        long dayInMillis = hourInMillis * 24;
        long yearInMillis = dayInMillis * 365;

        String elapsedYears = "0" + String.valueOf(diff / yearInMillis);
        diff = diff % yearInMillis;
        String elapsedDays = "0" + String.valueOf(diff / dayInMillis);
        diff = diff % dayInMillis;
        String elapsedHours = "0" + String.valueOf(diff / hourInMillis);
        diff = diff % hourInMillis;
        String elapsedMinutes = "0" + String.valueOf(diff / minuteInMillis);
        diff = diff % minuteInMillis;
        String elapsedSeconds = "0" + String.valueOf(diff / secondInMillis);

        String evaluationTime = elapsedHours.substring(elapsedHours.length()-2)  +":" + elapsedMinutes.substring(elapsedMinutes.length()-2) + ":" + elapsedSeconds.substring(elapsedSeconds.length()-2);
        
        System.out.println("Èas izvajanja = " + evaluationTime);
        population.printToFile(genFileName, "\nEvaluation time: " + evaluationTime);                

    }
}