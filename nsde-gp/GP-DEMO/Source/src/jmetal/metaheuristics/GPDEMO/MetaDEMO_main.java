//  MetaDEMO_main.java
//
//  Author:
//       Dejan Petelin <dejan.petelin@ijs.si>
//       Miha Mlakar <miha.mlakar@ijs.si>
//
//  Copyright (c) 2012 Dejan Petelin, Miha Mlakar
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
package jmetal.metaheuristics.GPDEMO;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

import com.mathworks.toolbox.javabuilder.MWException;

import jmetal.core.Algorithm;
import jmetal.core.Model;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.models.GP;
import jmetal.models.ModelFactory;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.problems.Kursawe;
import jmetal.problems.ProblemFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

/**
 * This class is the main program used to configure and run DEMO, a 
 * multiobjective scatter search metaheuristics, which is described in:
 *   ...
 */
public class MetaDEMO_main {
  public static Logger      logger_ ;      // Logger object
  public static FileHandler fileHandler_ ; // FileHandler object
    
  /**
  * @param args Command line arguments.
   * @throws Exception 
  */
  public static void main(String [] args) throws Exception {
    Problem   problem   ;         // The problem to solve
    Algorithm algorithm ;         // The algorithm to use
    Operator  selection ;
    Operator  crossover ;
    Model     model     ;
        
    HashMap  parameters ; // Operator parameters
    
    QualityIndicator indicators ; // Object to get quality indicators
        
    // Logger object and file to store log messages
    logger_      = Configuration.logger_ ;
    fileHandler_ = new FileHandler("Demo_main.log"); 
    logger_.addHandler(fileHandler_) ;
        
    model = null;
    problem = null;
    indicators = null;
    // problem
    if (args.length >= 1)
    	problem = (new ProblemFactory()).getProblem(args[0], new Object[]{"Real"});
    else
      problem = new Kursawe("Real", 3); 
    // model
    if (args.length >= 2)
    	model = (new ModelFactory()).getModel(args[1], problem.getNumberOfVariables(), problem.getNumberOfObjectives());
    else
      model = new GP(problem.getNumberOfVariables(), problem.getNumberOfObjectives());
    // indicators
    if (args.length >= 3)
    	indicators = new QualityIndicator(problem, args[2]);
        
    algorithm = new GPDEMO(problem, model);
    
    // Algorithm parameters
    algorithm.setInputParameter("populationSize", 100);
    algorithm.setInputParameter("maxIterations", 250);
    algorithm.setInputParameter("environmentalSelection", "s");
    
    // Crossover operator 
    parameters = new HashMap() ;
    parameters.put("CR", 0.5) ;
    parameters.put("F", 0.5) ;
    crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);   
    
    // Add the operators to the algorithm
    parameters = null ;
    selection = SelectionFactory.getSelectionOperator("DifferentialEvolutionSelection", parameters) ;

    algorithm.addOperator("crossover", crossover);
    algorithm.addOperator("selection", selection);
    
    // Execute the Algorithm 
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    long estimatedTime = System.currentTimeMillis() - initTime;
    
    // Result messages 
    logger_.info("Total execution time: " + estimatedTime + "ms");
    logger_.info("Variables values have been writen to file VAR");
    population.printVariablesToFile("VAR");    
    logger_.info("Objectives values have been writen to file FUN");
    population.printObjectivesToFile("FUN");

    if (indicators != null) {
      logger_.info("Quality indicators") ;
      logger_.info("Hypervolume: " + indicators.getHypervolume(population)) ;
      logger_.info("GD         : " + indicators.getGD(population)) ;
      logger_.info("IGD        : " + indicators.getIGD(population)) ;
      logger_.info("Spread     : " + indicators.getSpread(population)) ;
      logger_.info("Epsilon    : " + indicators.getEpsilon(population)) ;  
    } // if        
  } //main 
} //MetaDEMO_main
