//  ModelFactory.java
//
//  Author:
//       Dejan Petelin <dejan.petelin@ijs.si>
//
//  Copyright (c) 2011 Dejan Petelin
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

package jmetal.models;

import java.lang.reflect.Constructor;
import java.util.Properties;

import jmetal.core.Model;
import jmetal.core.Problem;
import jmetal.util.Configuration;
import jmetal.util.JMException;

/**
 * This class represents a factory for models
 */
public class ModelFactory {

	/**
   * Creates an object representing a model
   * @param name Name of the model
   * @return The object representing the model
   * @throws JMException 
   */
  public Model getModel(String name, int noInput, int noOutput) throws JMException {
    String base = "jmetal.models.";
    Model model;
    
    if (name == null)
    	return null;
    if (name.equalsIgnoreCase("SPGP"))
      name = "SPGPm";
    try {
      Class problemClass = Class.forName(base+name);
      Constructor constructors = problemClass.getConstructor(Properties.class);

      if (name.equalsIgnoreCase("SPGPm"))
      	model = (Model)constructors.newInstance(noInput, noOutput, 300);
      else
      	model = (Model)constructors.newInstance(noInput, noOutput);
      // TODO: Problem problem = (Model)constructors.newInstance(noInput, noOutput, params);
     
      return model;
    } // try
    catch(Exception e) {
      e.printStackTrace();
      Configuration.logger_.severe("ModelFactory.getModel: " +
        "Model '"+ name + "' does not exist. "  +
        "Please, check the model names in jmetal/models") ;
      throw new JMException("Exception in " + name + ".getModel()") ;
    } // catch
  }

}