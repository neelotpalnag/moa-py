//  SPEA2.java
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

package jmetal.metaheuristics.spea2;

import jmetal.util.JMException;
import jmetal.util.Ranking;
import jmetal.core.*;

import java.util.Comparator;
import java.util.Properties;
import jmetal.util.*;

/** 
 * This class representing the SPEA2 algorithm
 */
public class SPEA2 extends Algorithm{
          
  /**
   * Defines the number of tournaments for creating the mating pool
   */
  public static final int TOURNAMENTS_ROUNDS = 1;

  /**
  * Constructor.
  * Create a new SPEA2 instance
  * @param problem Problem to solve
  */
  public SPEA2(Problem problem) {                
    super(problem) ;
  } // Spea2
  
    /**
  * Constructor.
  * Create a new SPEA2 instance
  * @param problem Problem to solve
  */
  public SPEA2(Problem problem, Properties properties) {                
    super(problem, properties) ;
  } // Spea2
   
  /**   
  * Runs of the Spea2 algorithm.
  * @return a <code>SolutionSet</code> that is a set of non dominated solutions
  * as a result of the algorithm execution  
   * @throws JMException 
  */  
  public static SolutionSet environmentalSelection(SolutionSet union, int archiveSize){
      
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
  
  public SolutionSet execute() throws JMException, ClassNotFoundException {   
    int populationSize, archiveSize, maxEvaluations, evaluations;
    Operator crossoverOperator, mutationOperator, selectionOperator;
    SolutionSet solutionSet, archive, offSpringSolutionSet; 
    int frontGen;
    int frontFormat;
    
    String logFileName, frontFileName, genFileName, contents;
    
    //Read the params
    populationSize = Integer.parseInt(getParameter("Algorithm.populationSize"));
    archiveSize    = populationSize;
    maxEvaluations = Integer.parseInt(getParameter("Algorithm.maxEvaluations"));
    
    frontGen=Integer.parseInt(getParameter("Output.frontGen"));
    frontFormat=Integer.parseInt(getParameter("Output.frontFormat"));
    
    frontFileName=getParameter("Output.frontFileName");
    logFileName = getParameter("Output.logFileName");
    genFileName = getParameter("Output.genFileName");
        
    //Read the operators
    crossoverOperator = operators_.get("crossover");
    mutationOperator  = operators_.get("mutation");
    selectionOperator = operators_.get("selection");        
        
    //Initialize the variables
    // Pripravi prazno začetno populacijo (1.)
    solutionSet  = new SolutionSet(populationSize);
    // Pripravi prazen začetni arhiv (2.)
    archive     = new SolutionSet(archiveSize);
    // Postavi evaluations na 0 (3.)
    evaluations = 0;
        
    //-> Create the initial solutionSet
    // Naključno generiraj in ovrednosti zaćčetno populacijo(1.)
    Solution newSolution;
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);            
      problem_.evaluateConstraints(newSolution);
      evaluations++;
      solutionSet.add(newSolution);
    }                        
        
    // PONAVLJAJ (4.)
    while (evaluations < maxEvaluations){        
      
      // V nov arhiv prepiši vse nedominirane osebke iz P(t-1) in A(t-1) (4.2) 
      SolutionSet union = ((SolutionSet)solutionSet).union(archive);
      
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------    
      archive=environmentalSelection(union,archiveSize);         
      //-------------------------------------ENVIRONMENTAL SELECTION--------------------------------------------
      
      /*
       * Sicer populacijo P(t) generiraj iz arhiva A(t) z uporabo turnirske selekcije, križanja in mutacije (4.6)
       */                     
      // Create a new offspringPopulation
      offSpringSolutionSet= new SolutionSet(populationSize);    
      Solution  [] parents = new Solution[2];
      while (offSpringSolutionSet.size() < populationSize){           
        int j = 0;
        do{
          j++;                
          parents[0] = (Solution)selectionOperator.execute(archive);
        } while (j < SPEA2.TOURNAMENTS_ROUNDS); // do-while                    
        int k = 0;
        do{
          k++;                
          parents[1] = (Solution)selectionOperator.execute(archive);
        } while (k < SPEA2.TOURNAMENTS_ROUNDS); // do-while
            
        //make the crossover 
        Solution [] offSpring = (Solution [])crossoverOperator.execute(parents);
        //meke the mutation
        mutationOperator.execute(offSpring[0]);            
        problem_.evaluate(offSpring[0]);
        problem_.evaluateConstraints(offSpring[0]);            
        offSpringSolutionSet.add(offSpring[0]);
        evaluations++;
      } // while
      // End Create a offSpring solutionSet
      solutionSet = offSpringSolutionSet;
      
     /**
     * Write to the output file every frontGen generations
     */
      if(evaluations%frontGen == 0){
          solutionSet.printGenerationsToFile(frontFileName, frontFormat);      
      }
      solutionSet.printVariablesAndObjectivesToFile(logFileName);  
      
      contents=String.format("%d\t\t", evaluations);
      solutionSet.printToFile(genFileName, contents);
    } // while
    
    /*
     * Ovrednosti vse osebke iz populacije P(t) in arhiva A(t) (4.7)
     */
    Ranking ranking = new Ranking(archive);
    
    return ranking.getSubfront(0);
  } // execute    
} // SPEA2
