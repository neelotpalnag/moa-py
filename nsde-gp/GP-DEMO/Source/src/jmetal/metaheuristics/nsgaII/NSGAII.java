//  NSGAII.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
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

package jmetal.metaheuristics.nsgaII;

import java.util.Properties;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.core.*;
import jmetal.util.comparators.CrowdingComparator;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.*;

/**
 * This class implements the NSGA-II algorithm. 
 */
public class NSGAII extends Algorithm {

  /**
   * Constructor
   * @param problem Problem to solve
   */
  public NSGAII(Problem problem) {
    super (problem) ;
  } // NSGAII
  
  /**
   * Constructor
   * @param problem Problem to solve
   * @param properties Properties of algorithm
   */
  public NSGAII(Problem problem, Properties properties){
    super (problem, properties);  
  }

  /**   
   * Runs the NSGA-II algorithm.
   * @return a <code>SolutionSet</code> that is a set of non dominated solutions
   * as a result of the algorithm execution
   * @throws JMException 
   */
  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int populationSize;
    int maxEvaluations;
    int evaluations;
    int frontGen;
    int frontFormat;

    String frontFileName, logFileName, genFileName, contents;
    QualityIndicator indicators;    // QualityIndicator object
    int requiredEvaluations;        // Use in the example of use of the
                                    // indicators object (see below)
    SolutionSet population;
    SolutionSet offspringPopulation;
    SolutionSet union;

    Operator mutationOperator;
    Operator crossoverOperator;
    Operator selectionOperator;

    Distance distance = new Distance();

    //Read the parameters
    populationSize = Integer.parseInt(getParameter("Algorithm.populationSize"));
    maxEvaluations = Integer.parseInt(getParameter("Algorithm.maxEvaluations"));
    
    indicators = (QualityIndicator) getInputParameter("indicators");
    
    frontGen=Integer.parseInt(getParameter("Output.frontGen"));
    frontFormat=Integer.parseInt(getParameter("Output.frontFormat"));
    
    frontFileName=getParameter("Output.frontFileName");
    logFileName = getParameter("Output.logFileName");
    genFileName = getParameter("Output.genFileName");
    
    //Initialize the variables
    // 1. Generiranje populacije staršev P(0)
    population = new SolutionSet(populationSize);   
    evaluations = 0;

    requiredEvaluations = 0;

    //Read the operators
    mutationOperator = operators_.get("mutation");
    crossoverOperator = operators_.get("crossover");
    selectionOperator = operators_.get("selection");

    // Create the initial solutionSet
    // 1. Ovrednostimo začetno populacijo staršev
    Solution newSolution;
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      evaluations++;
      population.add(newSolution);
    } //for       
    
    population.printVariablesAndObjectivesToFile(logFileName);
    
    // Generations 
    while (evaluations < maxEvaluations) {  //Dokler ustavitveni kriterij ni izpolnjen
        
        
      // Create the offSpring solutionSet      
      // 2. Pripravimo prazno začetno populacijo potomcev Q(0)  
      offspringPopulation = new SolutionSet(populationSize);  
      
      // 4.8 Populacijo potomcev Q(t) generiraj iz populacije staršev P(t) z uporabo turnirske selekcije, križanja in mutacije
      Solution[] parents = new Solution[2];
      for (int i = 0; i < (populationSize / 2); i++) {
        if (evaluations < maxEvaluations) {
          //obtain parents
          parents[0] = (Solution) selectionOperator.execute(population);
          parents[1] = (Solution) selectionOperator.execute(population);
          Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
          mutationOperator.execute(offSpring[0]);
          mutationOperator.execute(offSpring[1]);
          problem_.evaluate(offSpring[0]);
          problem_.evaluateConstraints(offSpring[0]);
          problem_.evaluate(offSpring[1]);
          problem_.evaluateConstraints(offSpring[1]);
          offspringPopulation.add(offSpring[0]);
          offspringPopulation.add(offSpring[1]);
          evaluations += 2;
        } // if                            
      } // for

      // Create the solutionSet union of solutionSet and offSpring
      // 4.2 Združi stari populaciji staršev in potomcev
      union = ((SolutionSet) population).union(offspringPopulation);

      // Ranking the union 
      // 4.3 NE-DOMINIRANO SORTIRANJE
      Ranking ranking = new Ranking(union);

      int remain = populationSize;
      int index = 0;
      SolutionSet front = null;
      // 4.4 Pripravimo novo prazno populacijo staršev 
      population.clear();

      // Obtain the next front
      front = ranking.getSubfront(index);
      
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------
      // 4.5 V populacijo P daj prvih i front, ki še gredo cele v njo
      while ((remain > 0) && (remain >= front.size())) {
        //Assign crowding distance to individuals
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        //Add the individuals of this front
        for (int k = 0; k < front.size(); k++) {
          population.add(front.get(k));
        } // for

        //Decrement remain
        remain = remain - front.size();

        //Obtain the next front
        index++;
        if (remain > 0) {
          front = ranking.getSubfront(index);
        } // if        
      } // while

      // Remain is less than front(index).size, insert only the best one
      // 4.6 Fronto F(i+1), ki ne gre več cela v populacijo P(t) sortiraj z uporabo metrike nakopičenosti
      // 4.7 Populacijo P(t) dopolni z osebki iz F(i+1), ki so najmanj nakopičeni
      
      if (remain > 0) {  // front contains individuals to insert                        
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        front.sort(new CrowdingComparator());
        for (int k = 0; k < remain; k++) {
          population.add(front.get(k));
        } // for

        remain = 0;
      } // if                               
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------
            
      // This piece of code shows how to use the indicator object into the code
      // of NSGA-II. In particular, it finds the number of evaluations required
      // by the algorithm to obtain a Pareto front with a hypervolume higher
      // than the hypervolume of the true Pareto front.
      if ((indicators != null) && (requiredEvaluations == 0)) {
        double HV = indicators.getHypervolume(population);
        if (HV >= (0.98 * indicators.getTrueParetoFrontHypervolume())) {
          requiredEvaluations = evaluations;
        } // if
      } // if
    
      /**
       * Write to the output file every frontGen generations
       */
      if(evaluations%frontGen == 0){ 
    	//Mihas
         ranking.getSubfront(0).printVariablesAndObjectivesToFile(frontFileName);
         //population.printGenerationsToFile(frontFileName, frontFormat);       
      }
      
      population.printVariablesAndObjectivesToFile(logFileName);  
      
      contents=String.format("%d\t\t", evaluations);
      population.printToFile(genFileName, contents);

      
      
    } // while

    // Return as output parameter the required evaluations
    setOutputParameter("evaluations", requiredEvaluations);

    // Return the first non-dominated front
    // 4.9 Ovrednosti osebke iz populacije potomcev 
    Ranking ranking = new Ranking(population);
   
    return ranking.getSubfront(0);
  } // execute
} // NSGA-II
