//  Model.java
//
//  Authors:
//       Dejan Petelin <dejan.petelin@ijs.si>
//// 
//  Copyright (c) 2012 Dejan Petelin
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

package jmetal.core;

import java.io.Serializable;
import jmetal.util.JMException;

/** 
 *  This class implements a generic template for the model developed in
 *  jMetal. Every model must have ... 
 *  The class declares an abstract method called <code>execute</code>, which 
 *  defines the behavior of the algorithm.
 */ 
public abstract class Model implements Serializable {

  /**
   * Stores the number of inputs of the problem
   */
  protected int numberOfInputs_;
  
  /** 
   * Stores the number of outputs of the problem
   */
  protected int numberOfOutputs_;

  /**
   * Stores the model name
   */
  protected String modelName_ ;
  
  /**
   * Constructor
   */
  public Model(int noInput, int noOutput) {
  	numberOfInputs_ = noInput;
  	numberOfOutputs_ = noOutput;
  }

  /** 
   * Gets the number of decision inputs.
   * @return the number of inputs.
   */
  public int getNumberOfInputs() {
    return numberOfInputs_ ;   
  } // getNumberOfInputs
      
  /** 
   * Gets the the number of outputs.
   * @return the number of outputs.
   */
  public int getNumberOfOutputs() {
    return numberOfOutputs_ ;
  } // getNumberOfOutputs

  /**
   * Returns the problem name
   * @return The problem name
   */
  public String getName() {
    return modelName_ ;
  }

  /**
   * Returns the problem name
   * @return The problem name
   */
  public void setName(String name) {
    modelName_ = name ;
  }
  
  /**
   * Launches the execution of an specific algorithm.
   */
  public abstract void update(SolutionSet solutions) throws JMException, ClassNotFoundException;

  /**
   * Launches the execution of an specific algorithm.
   */
  public abstract void update(Solution solution) throws JMException, ClassNotFoundException;   

  /**
   * Launches the execution of an specific algorithm.
   */
  public abstract void evaluate(SolutionSet solutions) throws JMException, ClassNotFoundException;   

  /**
   * Launches the execution of an specific algorithm.
   */
  public abstract void evaluate(Solution solution) throws JMException, ClassNotFoundException;

  /**
   * Closes the connection to R for metrods that run through R.
   */
  public void closeRConnection()
  {
  	
  }

} // Model
