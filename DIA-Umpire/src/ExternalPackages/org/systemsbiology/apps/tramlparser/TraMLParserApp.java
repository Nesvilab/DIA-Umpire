/**
 * Copyright (c) 2010 Institute for Systems Biology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Contact info:
 * Mi-Youn K. Brusniak
 * Insitute for Systems Biology
 * 1441 North 34th St.
 * Seattle, WA  98103  USA
 * mbrusniak@systemsbiology.org
 *
 */
/**
 * Software: TraML Parser App version 0.1
 * Date: February 12, 2009
 * Very simple class that demonstrates the use of the ATAQS TraML parser. 
 *
 * Copyright (C) 2009 
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * David Campbell
 * Insitute for Systems Biology
 * 1441 North 34th St.
 * Seattle, WA  98103  USA
 * dcampbell@systemsbiology.org
 *
 */
package ExternalPackages.org.systemsbiology.apps.tramlparser;


// org.apache

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.LoggerContext;
import org.apache.logging.log4j.core.config.Configurator;
import org.apache.logging.log4j.core.config.DefaultConfiguration;

import java.nio.file.Paths;


/*  
 *
 */
public class TraMLParserApp {

    // Set a global logger and handler
    private static Logger logger = LogManager.getLogger(TraMLParserApp.class.getName());

    // Main class for cmd line invocation.
    public static void main(String[] args) {

        // MYKB: Use log4j properties
        String log4jPropertyFile = System.getProperty("log4j.configuration");
        if (log4jPropertyFile == null) {
            setDefaultLoggerHandler();
        } else {
//            PropertyConfigurator.configure(log4jPropertyFile);
            final LoggerContext context = (org.apache.logging.log4j.core.LoggerContext) LogManager.getContext(false);
            context.setConfigLocation(Paths.get(log4jPropertyFile).toUri());
        }

        TraMLParser parser = new TraMLParser();
        String transitionXML = "";
        // Interpret first arg as xml filename.
        try {
            parser.parse_file(args[0], logger);
            transitionXML = parser.getTransitionListXML(logger);
            parser.writeToFile("temp.txt");
        } catch (Exception e) {
            logger.error("Transport error: " + e.getMessage());
            e.printStackTrace();
        }
        // Get TraML document as a string.
        logger.info(transitionXML);
    } // End main

    private static void setDefaultLoggerHandler() {
        // Configure log4j style
        Configurator.initialize(new DefaultConfiguration());
        Configurator.setRootLevel(Level.INFO);

        // java.util.log configuration
        // logger.addHandler(new ConsoleHandler());

    } // End setDefaultLoggerHandler

} // End class
