//  DEMO.java
//
//  Author:
//       Miha Mlakar <miha.mlakar@ijs.si>
//       Dejan Ptelein <dejan.petelin@ijs.si>
//
//  Copyright (c) 2011 Miha Mlakar, Dejan Petelin
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
package jmetal.metaheuristics.demo;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.core.SolutionSet;
import jmetal.init.Params;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.util.Spea2Fitness;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.util.comparators.DominanceComparator;

/**
 * This class implements the DEMO algorithm. This algorithm is an adaptation
 * of the single-objective scatter search template defined by F. Glover in:
 * F. Glover. "A template for scatter search and path relinking", Lecture Notes 
 * in Computer Science, Springer Verlag, 1997. AbYSS is described in: 
 *   ...
 */
public class DEMO extends Algorithm{
  
  /**
  * Constructor
  * @param problem Problem to solve
  */
  public DEMO(Problem problem){
    super (problem);
  } // demo
  
  /**
  * Constructor
  * @param problem Problem to solve
  * @param properties Properties of algorithm
  */
  public DEMO(Problem problem, Properties properties){
    super (problem, properties);
  } // demo
  
  public SolutionSet environmentalSelectionSPEA2(SolutionSet union, int archiveSize){      
	  SolutionSet res;
	  /*
	   * Racunanje razdalj med osebki. Razdalje so shranjene v spremenljivki distance objekta Spa2Fitness
	  */
	  Spea2Fitness spea = new Spea2Fitness(union);
	  
	  /*
	   * Calculating:   strength(moc)
	   *                R(i) - rawFitness(groba uspesnost)
	   *                D(i) - distance (gostota)
	   *                F(i) = R(i) + D(i) - (prava uspesnost)
	   */
	  spea.fitnessAssign();
	  
	  /*
	   * Ce je nedominiranih osebkov vec kot je velikost novega arhiva, izloci osebke, ki so blizu drugim osebkom (4.3)
	   * Sicer pa, ce nedominirani osebki ne zapolnijo arhiva, arhiv dopolni z najboljsimi dominiranimi osebki iz P(t-1) in A(t-1) (4.4)
	   * Ce je zaustavitveni kritrij izpolnjen, koncaj (4.5)
	   */    
	  res = spea.environmentalSelection(archiveSize);
	  
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
  
  public void removeWorst(SolutionSet solutionSet, ArrayList<List<Double>> indicatorValues_, double maxIndicatorValue_) {
	  // Find the worst;
	  double worst = solutionSet.get(0).getFitness();
	  int worstIndex = 0;
	  double kappa = Double.parseDouble(getParameter("IBEA.kappa"));    //fitness scaling factor
	
	  for (int i = 1; i < solutionSet.size(); i++) {
      if (solutionSet.get(i).getFitness() > worst) {
	      worst = solutionSet.get(i).getFitness();
	      worstIndex = i;
      }
	  }
	
	  // Update the population
	  for (int i = 0; i < solutionSet.size(); i++) {
      if (i != worstIndex) {
	      double fitness = solutionSet.get(i).getFitness();
	      fitness -= Math.exp((-indicatorValues_.get(worstIndex).get(i) / maxIndicatorValue_) / kappa);
	      solutionSet.get(i).setFitness(fitness);
      }
	  }
	
	  // remove worst from the indicatorValues list
	  indicatorValues_.remove(worstIndex); // Remove its own list
	  Iterator<List<Double>> it = indicatorValues_.iterator();
	  while (it.hasNext()) {
      it.next().remove(worstIndex);
	  }
	
	  // remove the worst individual from the population
	  solutionSet.remove(worstIndex);
  } // removeWorst

  public Object[] calculateFitness(SolutionSet solutionSet) {
	  // Obtains the lower and upper bounds of the population
	  double[] maximumValues = new double[problem_.getNumberOfObjectives()];
	  double[] minimumValues = new double[problem_.getNumberOfObjectives()];
	
	  for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
      maximumValues[i] = -Double.MAX_VALUE; // i.e., the minus maxium value
      minimumValues[i] = Double.MAX_VALUE; // i.e., the maximum value
	  }
	
	  for (int pos = 0; pos < solutionSet.size(); pos++) {
      for (int obj = 0; obj < problem_.getNumberOfObjectives(); obj++) {
	      double value = solutionSet.get(pos).getObjective(obj);
	      if (value > maximumValues[obj]) {
	      	maximumValues[obj] = value;
	      }
	      if (value < minimumValues[obj]) {
	      	minimumValues[obj] = value;
	      }
      }
	  }
	
	  Object[] result = computeIndicatorValuesHD(solutionSet, maximumValues, minimumValues);
	  for (int pos = 0; pos < solutionSet.size(); pos++) {
      fitness(solutionSet, pos, (ArrayList<List<Double>>)result[0], Double.valueOf(result[1].toString()));
	  }
	  
	  return result;
  }

  /**
   * This structure store the indicator values of each pair of elements
   */
  public Object[] computeIndicatorValuesHD(SolutionSet solutionSet, double[] maximumValues, double[] minimumValues) {
	  SolutionSet A, B;
	  // Initialize the structures
	  ArrayList<List<Double>> indicatorValues_ = new ArrayList<List<Double>>();
	  double maxIndicatorValue_ = -Double.MAX_VALUE;
	
	  for (int j = 0; j < solutionSet.size(); j++) {
      A = new SolutionSet(1);
      A.add(solutionSet.get(j));

      List<Double> aux = new ArrayList<Double>();
      for (int i = 0; i < solutionSet.size(); i++) {
	      B = new SolutionSet(1);
	      B.add(solutionSet.get(i));
	
	      int flag = (new DominanceComparator()).compare(A.get(0), B.get(0));
	
	      double value = 0.0;
	      if (flag == -1) {
          value = -calcHypervolumeIndicator(A.get(0), B.get(0), problem_.getNumberOfObjectives(), maximumValues, minimumValues);
	      } else {
          value = calcHypervolumeIndicator(B.get(0), A.get(0), problem_.getNumberOfObjectives(), maximumValues, minimumValues);
	      }
	      //double value = epsilon.epsilon(matrixA,matrixB,problem_.getNumberOfObjectives());

	      //Update the max value of the indicator
	      if (Math.abs(value) > maxIndicatorValue_) {
          maxIndicatorValue_ = Math.abs(value);
	      }
	      aux.add(value);
      }
      indicatorValues_.add(aux);
	  }
	  return new Object[]{indicatorValues_, maxIndicatorValue_};
  } // computeIndicatorValues

  /**
   * Calculate the fitness for the individual at position pos
   */
  public void fitness(SolutionSet solutionSet, int pos, ArrayList<List<Double>> indicatorValues_, double maxIndicatorValue_) {
	  double fitness = 0.0;
	  double kappa = Double.parseDouble(getParameter("IBEA.kappa"));
	
	  for (int i = 0; i < solutionSet.size(); i++) {
      if (i != pos) {
      	fitness += Math.exp((-1 * indicatorValues_.get(i).get(pos) / maxIndicatorValue_) / kappa);
      }
	  }
	  solutionSet.get(pos).setFitness(fitness);
  }

  double calcHypervolumeIndicator(Solution p_ind_a, Solution p_ind_b, int d, double maximumValues[], double minimumValues[]) {
	  double a, b, r, max;
	  double volume = 0;
	  double rho = Double.parseDouble(getParameter("IBEA.rho")); 
	
	  r = rho * (maximumValues[d - 1] - minimumValues[d - 1]);
	  max = minimumValues[d - 1] + r;
	
	  a = p_ind_a.getObjective(d - 1);
	  if (p_ind_b == null) {
	    b = max;
	  } else {
	    b = p_ind_b.getObjective(d - 1);
	  }
	
	  if (d == 1) {
      if (a < b) {
      	volume = (b - a) / r;
      }
      else {
      	volume = 0;
      }
	  } else {
      if (a < b) {
			  volume = calcHypervolumeIndicator(p_ind_a, null, d - 1, maximumValues, minimumValues) * (b - a) / r;
			  volume += calcHypervolumeIndicator(p_ind_a, p_ind_b, d - 1, maximumValues, minimumValues) * (max - b) / r;
      }
      else {
      	volume = calcHypervolumeIndicator(p_ind_a, p_ind_b, d - 1, maximumValues, minimumValues) * (max - b) / r;
      }
	  }
	
	  return (volume);
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
        
                
        // Create the initial solutionSet
        Solution newSolution;
        for (int i = 0; i < populationSize; i++) {
          newSolution = new Solution(problem_);                    
          problem_.evaluate(newSolution);            
          problem_.evaluateConstraints(newSolution);          
          evaluations++;
          if(logMode!=0)
        	  newSolution.printVariablesandObjectivesToFile(logFileName, evaluations);
          population.add(newSolution);
        } //for       
       
        if(frontMode != 0) {
        ranking = new Ranking(population);
        if(frontFormat==1) //decision and objective vector
      	  ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
        else if(frontFormat==0) //only objective vector
      	  ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);  
        }
        
        //first iteration is initial population
        iterations ++;
        
        //printing into the generations file
        if(genMode!=0) //0 = no output
        	population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, -1, -1, -1, 0, -1.0);
                     
        
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

                    /*
                     * Repair the candidate if it falls out of bounds of the decision space
                     */
                    
                    
                    //Candidate evaluation
                    problem_.evaluate(child) ;
                    problem_.evaluateConstraints(child);
                    evaluations++ ;
                    if(logMode!=0)
                    	child.printVariablesandObjectivesToFile(logFileName, evaluations);
                 	                    
                    // Dominance test
                    int result  ;
                    if(areSolutionsTheSame(population.get(i), child) == false)
                    {
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
            	  while (population.size() > populationSize) {		//offspringPopulation
                      Object[] result = calculateFitness(population);
                      ArrayList<List<Double>> indicatorValues_ = (ArrayList<List<Double>> )result[0];
                      double maxIndicatorValue_ = Double.valueOf(result[1].toString());
                      removeWorst(population, indicatorValues_, maxIndicatorValue_);
                  }   
              }else if(selectionProcedure==2){   //SPEA 2
                  population=environmentalSelectionSPEA2(population, populationSize);	//offspringPopulation
              }              
                          
              /* we get population of reduced (half) size after the environmental selection               
               * so we double it back. 
               */
              population.setMaxSize(populationSize*2);
              
              
              
             /**
              * Write to the fronts output file every frontGen generations
              */
              if(iterations%frontGen == 0){           	     
            	  if(frontMode != 0) {
            		  population.printToFile(frontFileName, "\n");
            		  //Ranking of solutions	              
	                  ranking = new Ranking(population);
	                  if(frontFormat==1) //decision and objective vector
	                	  ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);   
	                  else if(frontFormat==0) //only the objective vector
	                	  ranking.getSubfront(0).printGenerationsToFile(frontFileName, 0);
            	  }
              }                       
              
              iterations ++ ;
             
              if(genMode!=0) //0 = no output 
            	  population.printInfoAboutGenerations(genFileName, problem_, population, populationSize, iterations, numberOfChildBetter, numberOfParentBetter, numberOfIncomparable, 0, -1.0);
              
            } // while

            // Return the first non-dominated front
            ranking = new Ranking(population);        
            return ranking.getSubfront(0);
       }// execute

 private boolean areSolutionsTheSame(Solution parent, Solution child) throws JMException {
  		
  		for(int i=0; i<Params.getNumberOfVariables(); i++)
  		{
  			if(parent.getDecisionVariables()[i].getValue() != child.getDecisionVariables()[i].getValue())
  				return false;
  		}
  		
  		return true;
  	}
    
} // demo