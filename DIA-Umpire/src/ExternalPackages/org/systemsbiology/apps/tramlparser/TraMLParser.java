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
 * Software: TraML Parser version 0.1
 * Date: February 12, 2009
 * Reader and writer of TraML schema using xmlbean classes
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

// Core

import org.apache.logging.log4j.Logger;
import ExternalPackages.org.hupo.psi.ms.traml.TraMLType;
import ExternalPackages.org.systemsbiology.constants.JTRAML_URL;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;
import javax.xml.namespace.QName;
import java.io.File;
import java.io.FileWriter;
import java.io.StringWriter;

// third-party

public class TraMLParser {

    // Define instance variables.
    private TraMLType traML;


    // Constructor
    public TraMLParser() {

    }

    /**
     * Set TraML object for parsing
     *
     * @param newTraML
     */
    public void setTraML(TraMLType newTraML) {
        this.traML = newTraML;
    }

    /**
     * Get TraML object
     */
    public TraMLType getTraML() {
        return this.traML;
    }

    /**
     * Get parsed TraML object
     *
     * @param tramlFileName - file name
     * @param logger        - logger
     * @throws Exception - exception during parsing
     */
    public void parse_file(String tramlFileName, Logger logger) throws Exception {

        JAXBContext ctx = JAXBContext.newInstance("ExternalPackages.org.hupo.psi.ms.traml");
        Unmarshaller um = ctx.createUnmarshaller();
        JAXBElement<TraMLType> jaxb_traml = (JAXBElement<TraMLType>) um.unmarshal(new File(tramlFileName));
        this.traML = jaxb_traml.getValue();

    } // End parse_file

    public String getTransitionListXML(Logger logger) throws Exception {

        JAXBContext ctx = JAXBContext.newInstance("ExternalPackages.org.hupo.psi.ms.traml");
        Marshaller m = ctx.createMarshaller();
        m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        m.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, JTRAML_URL.TRAML_XSD_LOCATION);
        m.setProperty(Marshaller.JAXB_ENCODING, "UTF-8");

        JAXBElement<TraMLType> tramlWrap =
                new JAXBElement<TraMLType>(new QName(JTRAML_URL.TRAML_URI, "TraML"), TraMLType.class, traML);

        StringWriter sw = new StringWriter();
        m.marshal(tramlWrap, sw);

        return sw.toString();
    } // End getTransitionListXML

    public void writeToFile(String filename) throws Exception {

        JAXBContext ctx = JAXBContext.newInstance("ExternalPackages.org.hupo.psi.ms.traml");
        Marshaller m = ctx.createMarshaller();
        m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
        m.setProperty(Marshaller.JAXB_SCHEMA_LOCATION, JTRAML_URL.TRAML_XSD_LOCATION);
        m.setProperty(Marshaller.JAXB_ENCODING, "UTF-8");

        JAXBElement<TraMLType> tramlWrap =
                new JAXBElement<TraMLType>(new QName(JTRAML_URL.TRAML_URI, "TraML"), TraMLType.class, traML);

        m.marshal(tramlWrap, new FileWriter(filename));
    }
} // End class

