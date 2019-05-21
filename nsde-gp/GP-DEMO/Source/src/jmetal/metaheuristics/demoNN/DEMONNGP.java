/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package jmetal.metaheuristics.demoNN;

import java.io.FileReader;
import java.util.Comparator;
import java.util.Properties;
import java.util.Random;

import org.omg.CORBA.portable.IndirectionException;

import weka.associations.gsp.Element;
import weka.classifiers.Evaluation;
import weka.classifiers.functions.LibSVM;
import weka.classifiers.functions.LinearRegression;
import weka.classifiers.functions.MultilayerPerceptron;
import weka.classifiers.pmml.consumer.NeuralNetwork;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.pmml.MiningSchema;
import weka.filters.supervised.instance.Resample;

import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.init.Params;
import jmetal.metaheuristics.spea2.SPEA2;
import jmetal.qualityIndicator.Hypervolume;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.Spea2Fitness;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;

/**
 *
 * @author Rok
 */
public class DEMONNGP extends Algorithm{
    
	private SolutionSet setOfExactlyEvaluatedSolutions;
	private Model modelGlobal_;
	
	
 /**
  * Constructor
  * @param problem Problem to solve
  */
  public DEMONNGP(Problem problem){
    super (problem) ;
  } // demo
  
 /**
  * Constructor
  * @param problem Problem to solve
  */
  public DEMONNGP(Problem problem, Model globalModel, Model localModel, Properties properties){
    super (problem, properties);

   	modelGlobal_ = globalModel;
    
  } // demo
  
  public SolutionSet environmentalSelectionSPEA2(SolutionSet union, int archiveSize){
      
      SolutionSet res;
      /*
       * Računanje razdalj med osebki. Razdalje so shranjene v spremenljivki distance objekta Spa2Fitness
      */
      Spea2Fitness spea = new Spea2Fitness(union);
      
      /*
       * Calculating:   strength(moč)
       *                R(i) - rawFitness(groba uspešnost)
       *                D(i) - distance (gostota)
       *                F(i) = R(i) + D(i) - (prava uspešnost)
       */
      spea.fitnessAssign();
      
      /*
       * Če je nedominiranih osebkov več kot je velikost novega arhiva, izloči osebke, ki so blizu drugim osebkom (4.3)
       * Sicer pa, če nedominirani osebki ne zapolnijo arhiva, arhiv dopolni z najboljšimi dominiranimi osebki iz P(t-1) in A(t-1) (4.4)
       * Če je zaustavitveni kritrij izpolnjen, končaj (4.5)
       */    
      res=spea.environmentalSelection(archiveSize);  
      
      return res;
  }

  public SolutionSet environmentalSelectionNSGAII(Ranking ranking, Distance distance, SolutionSet front, int remain, int index){
      
        SolutionSet finalPop=new SolutionSet(remain);

        // Obtain the next front
        front=ranking.getSubfront(index);

        while ((remain > 0) && (remain >= front.size())){ 
            //Assign crowding distance to individuals
            distance.crowdingDistanceAssignment(front,problem_.getNumberOfObjectives());

            for (int k = 0; k < front.size(); k++ ) {
                finalPop.add(front.get(k));
            } // for

            //Decrement remain
            remain = remain - front.size();

            //Obtain the next front
            index++;
            if (remain > 0) {
                front = ranking.getSubfront(index);
            } // if
         }//end while

         // remain is less than front(index).size, insert only the best one
         if (remain > 0) {  // front contains individuals to insert                        
            while (front.size() > remain) {
                distance.crowdingDistanceAssignment(front,problem_.getNumberOfObjectives());
                front.remove(front.indexWorst(new CrowdingComparator()));
            }
            for (int k = 0; k < front.size(); k++) {
                finalPop.add(front.get(k));
            }
            remain = 0; 
          } // if    
      
      
      return finalPop;
      
  }
  
    @Override
  public SolutionSet execute() throws JMException, ClassNotFoundException {
        int populationSize ;
        int maxEvaluations  ;
        int maxIterations;
        int evaluations    ;
        int iterations     ;
        int selectionProcedure;
        int frontGen;
        int frontFormat;
        int frontMode;
        int logMode;
        int genMode;
        Ranking ranking;
        int numberOfFeasible;
        int numberOfInfeasible;
        int numberOfNonCalculated;
        double constraintViolation;
        int numberOfChildBetter;
        int numberOfParentBetter;
        int numberOfIncomparable;

        String logFileName, frontFileName, genFileName, contents;
        
        SolutionSet population          ;

        Distance   distance  ;
        Comparator dominance ;
        
        Operator selectionOperator ;
        Operator crossoverOperator ;
        
        distance  = new Distance()  ;               
        dominance = new DominanceComparator();
        
        Solution parents[] ;
        
        //Read the parameters
        //Mihas
        populationSize = Integer.parseInt(getParameter("Algorithm.populationSize"));
        selectionProcedure = Integer.parseInt(getParameter("Algorithm.selectionProcedure"));
        maxEvaluations = Integer.parseInt(getParameter("Algorithm.maxEvaluations"));
        maxIterations = maxEvaluations/populationSize + (int)Math.ceil(maxEvaluations%populationSize);
        //populationSize = Integer.parseInt(getParameter("populationSize"));
        //maxIterations  = Integer.parseInt(getParameter("maxIterations"));
        
        frontGen=Integer.parseInt(getParameter("Output.frontGen"));
        frontFormat=Integer.parseInt(getParameter("Output.frontFormat"));
        frontMode = Integer.parseInt(getParameter("Output.frontMode"));
    
        frontFileName=getParameter("Output.frontFileName");
        logFileName = getParameter("Output.logFileName");
        genFileName = getParameter("Output.genFileName");
        genMode=Integer.parseInt(getParameter("Output.genMode"));
        logMode = Integer.parseInt(getParameter("Output.logMode"));
        
        selectionOperator = operators_.get("selection");   
        crossoverOperator = operators_.get("crossover");
        
        //Initialize the variables
        population  = new SolutionSet(populationSize * 2);        
        evaluations = 0;                
        iterations  = 0;      
        final int numberOfExactEvaluationsOfGenetaions = Params.getNumberOfExactEvaluations();
        int currentGenerationOfExactEvaluations = 0;
        final int numberOfApproximateEvaluationsOfGenetaions = Params.getNumberOfApproximateEvaluations();
        int currentGenerationOfApproximateEvaluations = 0;
        boolean useExactEvaluation = true;
        int additionalEvaluations =0;
        double NNAccuracy = -1;
        double epsilon = 0.1;

        setOfExactlyEvaluatedSolutions = new SolutionSet(Integer.MAX_VALUE);
        
        
        // Create the initial solutionSet
        Solution newSolution;
        for (int i = 0; i < populationSize; i++) {
          newSolution = new Solution(problem_);                    
          problem_.evaluate(newSolution);            
          problem_.evaluateConstraints(newSolution); 
          newSolution.setExactllyEvaluated(true);
          setOfExactlyEvaluatedSolutions.add(newSolution);
          evaluations++;
          if(logMode!=0)
        	  newSolution.printVariablesandObjectivesToFile(logFileName, evaluations);
          population.add(newSolution);
        } //for       
       
        //first iteration is initial population
        iterations ++;
        currentGenerationOfExactEvaluations ++; 	//first generation is always exact
        if(currentGenerationOfExactEvaluations >= numberOfExactEvaluationsOfGenetaions)
        {
        	currentGenerationOfExactEvaluations = 0;
        	useExactEvaluation = false;
        	
        	modelGlobal_.update(setOfExactlyEvaluatedSolutions);
        	setOfExactlyEvaluatedSolutions.clear();
        }
        
        //print front
        if(iterations%frontGen == 0){   
	        if(frontMode != 0) {
	        	ranking = new Ranking(population);
	        	if(frontFormat==1) //decision and objective vector
	        		ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
	        	else if(frontFormat==0) //only objective vector
	        		ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);  
	        }
        }
  
        //printing into the generations file
        if(genMode!=0) //0 = no output
        	if(currentGenerationOfExactEvaluations == 0 && (currentGenerationOfApproximateEvaluations >0  || useExactEvaluation == true) )
      		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, -1, -1, -1, additionalEvaluations, NNAccuracy);
      	  else
      		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, -1, -1, -1, additionalEvaluations, -1);	      	
          
        
        
        // Generations ...
        while (iterations < maxIterations) {
              // Create the offSpring solutionSet (empty)      
              //offspringPopulation  = new SolutionSet(populationSize*2);        

        	    numberOfChildBetter = 0;
              numberOfParentBetter = 0;
              numberOfIncomparable = 0;	
              for (int i = 0; i < populationSize; i++){   
                    // Obtain parents. Two parameters are required: the population and the 
                    //                 index of the current individual
                    parents = (Solution [])selectionOperator.execute(new Object[] { population, i} );            //GENERIRANJE KANDIDATA

                    Solution child ;
                    // Crossover. Two parameters are required: the current individual and the 
                    //            array of parents
                    child = (Solution)crossoverOperator.execute(new Object[]{population.get(i), parents});   
                    
                    //Candidate evaluation
                    if(useExactEvaluation)
                    {
                    	problem_.evaluate(child) ;
                    	problem_.evaluateConstraints(child);
                    	child.setExactllyEvaluated(true);
                    	setOfExactlyEvaluatedSolutions.add(child);
                    }
                    else
                    {
                    	approximateSolution(child);                    	
                    	problem_.evaluateConstraints(child);
                    	child.setExactllyEvaluated(false);
                    }
                    
                    evaluations++ ;
                    if(logMode!=0)
                    	child.printVariablesandObjectivesToFile(logFileName, evaluations);
                 	                    
                    // Dominance test
                    if(areSolutionsTheSame(population.get(i), child) == false)
	                    {
	                    int result  ;
	                    result = dominance.compare(population.get(i), child) ;
	                    if (result == -1) { // Solution i dominates child
	                    	numberOfParentBetter ++;                    
	                    } // if
	                    else if (result == 1) { // child dominates
	                    	population.replace(i, child);
	                    	numberOfChildBetter ++;                    	                      
	                    } // else if
	                    else { // the two solutions are non-dominated
	                    	population.add(child) ;
	                        numberOfIncomparable ++;
	                    } // else
                    }
                    else // the two solutions are the same so we do not add the child
                    {
                    	numberOfIncomparable ++;
                    }
              } // for           

          
              /* Environmental selection
               * Demo modifikacija: PARETO OPTIMALNE FRONTE + METRIKE NAKOPIČENOSTI
               * Demo mora omogočat, da to delamo na 3 načine in sicer: - NSGA 2
               *                                                        - IBEA
               *                                                        - SPEA 2                                                        
               */
              
              if(selectionProcedure==0){         //NSGA 2
            	// Ranking the offspring population
                  ranking = new Ranking(population); 
                  int remain = populationSize;
                  int index  = 0;
                  SolutionSet front = null;
                  population=environmentalSelectionNSGAII(ranking, distance, front, remain, index);
              }
              else if(selectionProcedure==1){    //IBEA
                  while(population.size() > populationSize){		//offspringPopulation
                      //removeWorst(offspringPopulation);
                  }
                  // Next    
              }else if(selectionProcedure==2){   //SPEA 2
                  population=environmentalSelectionSPEA2(population, populationSize);	//offspringPopulation
              }              
                          
              /* we get population of reduced (half) size after the environmental selection               
               * so we double it back. 
               */
              population.setMaxSize(populationSize*2);                                                             
              
              iterations ++ ;
              //decide wheather to use exact or approximate evaluation in the next generation
              if(useExactEvaluation)
              {
            	  currentGenerationOfExactEvaluations ++;
            	  if(currentGenerationOfExactEvaluations >= numberOfExactEvaluationsOfGenetaions)
                  {
                  	 currentGenerationOfExactEvaluations = 0;
                  	 useExactEvaluation = false;
                  	 modelGlobal_.update(setOfExactlyEvaluatedSolutions);
                  }
              }
              else	//use meta models
              {
            	  currentGenerationOfApproximateEvaluations ++;
            	  if(currentGenerationOfApproximateEvaluations >= numberOfApproximateEvaluationsOfGenetaions)
                {
            		  currentGenerationOfApproximateEvaluations = 0;
            		  useExactEvaluation = true;
            		  
            		  //reevaluate NN front solutions 
            		  int tempAdditionalEvaluations;
            		  
            		  if(1==2)
            		  {
	            		  do
	            		  {
	            			  tempAdditionalEvaluations =0;
	            			  ranking = new Ranking(population);
	            			  for(int i=0;i<ranking.getSubfront(0).size();i++)
	                		  {
	                			  if(ranking.getSubfront(0).get(i).isExactllyEvaluated() == false)   //this means that it was calculated by NN
	                			  {
	                				  additionalEvaluations++;
	                				  tempAdditionalEvaluations++;
	                				  problem_.evaluate(ranking.getSubfront(0).get(i));
	                				  problem_.evaluateConstraints(ranking.getSubfront(0).get(i));
	                				  ranking.getSubfront(0).get(i).setExactllyEvaluated(true);
	                			  }            			  
	                		  }  
	            			  
	            		  } while(tempAdditionalEvaluations > 0);            		              		              		  
            		  }
            		  //NNAccuracy = -1; 		//that we write this only when it was calculated
                }
              }
             
              
              /**
               * Write to the fronts output file every frontGen generations
               */
               if(iterations%frontGen == 0){           	     
             	  if(frontMode != 0) {
             		  population.printToFile(frontFileName, "\n");		//print new line
             		  //Ranking of solutions	              
 	                  ranking = new Ranking(population);
 	                  if(frontFormat==1) //decision and objective vector
 	                	  ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
 	                  else if(frontFormat==0) //only objective vector
 	                	  ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);
             	  }
               } 
               
              //print generations output
              if(genMode!=0) //0 = no output
              {
            	  if(currentGenerationOfExactEvaluations == 0 && (currentGenerationOfApproximateEvaluations >0  || useExactEvaluation == true) )
            		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, numberOfChildBetter, numberOfParentBetter, numberOfIncomparable, additionalEvaluations, NNAccuracy);
            	  else
            		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, numberOfChildBetter, numberOfParentBetter, numberOfIncomparable, additionalEvaluations, -1);
              }
            
            } // while

        
        //izra�umano exzaktno re�itve na fronti
        int tempAdditionalEvaluations = -1;
        int finalAdditionalEvaluations = 0;
        if(1==1)
  		  {
        	
    		  do
    		  {
    			  tempAdditionalEvaluations =0;
    			  ranking = new Ranking(population);
    			  for(int i=0;i<ranking.getSubfront(0).size();i++)
        		  {
        			  if(ranking.getSubfront(0).get(i).isExactllyEvaluated() == false)   //this means that it was calculated by NN
        			  {        				  
        			  	finalAdditionalEvaluations++;
        			  	tempAdditionalEvaluations++;
        				  problem_.evaluate(ranking.getSubfront(0).get(i));
        				  problem_.evaluateConstraints(ranking.getSubfront(0).get(i));
        				  ranking.getSubfront(0).get(i).setExactllyEvaluated(true);
        			  }            			  
        		  }  
    			  
    		  } while(tempAdditionalEvaluations > 0);
    		  
    		  System.out.println("Temp dodatni izracuni " + finalAdditionalEvaluations);
    		  ranking = new Ranking(population);  
    		  if(frontFormat==1) //decision and objective vector
         	  ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
           else if(frontFormat==0) //only objective vector
         	  ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);
    		  
    		  if(genMode!=0) //0 = no output
          {
        	  if(currentGenerationOfExactEvaluations == 0 && (currentGenerationOfApproximateEvaluations >0  || useExactEvaluation == true) )
        		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations + 1000, -1, -1, -1, finalAdditionalEvaluations, NNAccuracy);
        	  else
        		  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations + 1000, -1, -1, -1, finalAdditionalEvaluations, -1);
          }
  		  }
        
        else
        {
            // Return the first non-dominated front
        System.out.println("Dodatni izracuni " + additionalEvaluations);
            ranking = new Ranking(population);  
        }
        
            return ranking.getSubfront(0);

       }//end execute  
    
    private void approximateSolution(Solution child) throws ClassNotFoundException, JMException {

    	double distanceToFisibility;
  	
    	modelGlobal_.evaluate(child);
    	
    	child.setExactllyEvaluated(false);
    	
    	
    	if(problem_.getName() == "StoreSteel" || problem_.getName() == "ExternEvaluator")
    	{
    		
  	  	if(problem_.getNumberOfFeatures() > 0)
  	  	{	
  	  		//calculate objectives from features
  	  		for(int i=0; i<problem_.getNumberOfFeatures(); i++)
  	  			child.setObjective(i, getObjective(i, child.getFeature(i)));  		
  	  	}  	
    	
  	  	distanceToFisibility =  calculateDistanceToFeasibility(child);
  	  	if(distanceToFisibility != 0)
  	  	{
  	  		child.setNumberOfViolatedConstraint(1);
  	  		child.setOverallConstraintViolation(distanceToFisibility);
  	  	}
  	  	
  	  	else
  	  	{
  	  		child.setNumberOfViolatedConstraint(0);
  	  		child.setOverallConstraintViolation(0);
  	  	}
    	}
    }

    
    private double calculateDistanceToFeasibility(Solution solution) {
    	double distance = 0;
    			
    			if(solution.getFeature(0) < 10)
    				distance = distance + 10 - solution.getFeature(0);
    			else if(solution.getFeature(0) > 11)
    				distance = distance + solution.getFeature(0) - 11;
    			
    			if(solution.getFeature(1) < 11)
    				distance = distance + 11 - solution.getFeature(1);
    			else if(solution.getFeature(1) > 15)
    				distance = distance + (solution.getFeature(1) - 15)/4;
    			
    			if(solution.getFeature(2) < 1115)
    				distance = distance + 1115 - solution.getFeature(2);
    			else if(solution.getFeature(2) > 1130)
    				distance = distance + (solution.getFeature(2) - 1130)/15;
    			
    		  return (-1 * distance);
    	  }

    		private double getObjective(int i, double feature) {

    			double objective;
    		  
    		  if(i== 0)
    		  	objective = Math.abs(feature-10);
    		  else if(i== 1)
    		  	objective = Math.abs(feature-13);
    		  else
    		  	objective = Math.abs(feature-1122.5);
    		  
    		  return objective;
    	  }
    		
	private boolean areSolutionsTheSame(Solution parent, Solution child) throws JMException {
		
		for(int i=0; i<Params.getNumberOfVariables(); i++)
		{
			if(parent.getDecisionVariables()[i].getValue() != child.getDecisionVariables()[i].getValue())
				return false;
		}
		
		return true;
	}

	private Instances[] updateExamplesDataBase(Instances[] iExamples, Solution newSolution) {
    	//For every output value we need to create new neural network
      	for(int i=0; i<Params.getNumberOfObjectives(); i++)
      	{      		
      		try
            {           			
      			Instance iExample;      			      			
      			iExample = new Instance(Params.getNumberOfVariables() + 1);
      			for(int k=0; k<Params.getNumberOfVariables();k++)      				
      				iExample.setValue(k, newSolution.getDecisionVariables()[k].getValue());
      					
      			iExample.setValue(Params.getNumberOfVariables(), newSolution.getObjective(i));
      			       			 
      			// add the instance
      			iExamples[i].add(iExample);      			       
                    	            	
      		} 
            catch (Exception e) 
            { 
            	e.printStackTrace(); 	
            }       
      	}
		return iExamples;
	}

	private Instances getInitializedInstances() {

        FastVector attr = new FastVector(Params.getNumberOfVariables() + 1);	//+1 for the criterion
        Attribute Attribute1;
        for(int i=0; i< Params.getNumberOfVariables(); i++)
        {
        	Attribute1 = new Attribute(i+ "- input variable");
        	attr.addElement(Attribute1);
        }
        
        Attribute1 = new Attribute("- class");
    	attr.addElement(Attribute1);
    	
    	
		
		Instances iExamples = new Instances("Examples", attr, 0); 
		iExamples.setClassIndex(Params.getNumberOfVariables());
		return iExamples;
	}

	public static Instances getInstances(String[] attributes, String[] classes, double[][] data,
			String[] dataClasses) {

		// J48 tree;
		FastVector atts, attVals;
		Instances learnData;
		double[] vals;

		atts = new FastVector();
		for (int i = 0; i < attributes.length; i++) {
			atts.addElement(new weka.core.Attribute(attributes[i], weka.core.Attribute.NUMERIC));
		}
		attVals = new FastVector();
		for (int i = 0; i < classes.length; i++) {
			attVals.addElement(classes[i] + "");
		}

		atts.addElement(new weka.core.Attribute("state", attVals));

		learnData = new Instances("Learn_data", atts, 0);
		for (int i = 0; i < data.length; i++) {
			vals = new double[data[i].length + 1];
			for (int j = 0; j < data[i].length; j++) {
				vals[j] = data[i][j] + 0.0;
			}
			vals[data[i].length] = attVals.indexOf(dataClasses[i] + "");
			
			learnData.add(new Instance(1.0, vals));			
		}

		try {
			Resample resample = new Resample();
			String[] resampleOptions = new String[1];
			resampleOptions[0] = "-B 1";
			resample.setOptions(resampleOptions);
			resample.setInputFormat(learnData);
			learnData = Resample.useFilter(learnData, resample);
		}
		catch (Exception e) {
			System.out.println(" problem with resample " + e.toString());
		}

		learnData.setClassIndex(learnData.numAttributes() - 1);
		return learnData;
	}

}//end demo
