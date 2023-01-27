/*
 * This file is part of DIA-Umpire.
 *
 * DIA-Umpire is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * DIA-Umpire is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with DIA-Umpire. If not, see <https://www.gnu.org/licenses/>.
 *
 */
package MSUmpire.SpectrumParser;

import MSUmpire.BaseDataStructure.SpectralDataType;
import MSUmpire.BaseDataStructure.XYData;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.TreeMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Definitions of DIA type and isolation window settings
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class DIA_Setting implements Serializable{
    private static final long serialVersionUID = 646984181894L;
    public TreeMap<XYData, ArrayList<Integer>> DIAWindows = new TreeMap<XYData, ArrayList<Integer>>();
    public TreeMap<XYData, ArrayList<Integer>> MS1Windows = new TreeMap<XYData, ArrayList<Integer>>();
    public float F_DIA_WindowSize = 25;
    public SpectralDataType.DataType dataType;
 
    public void WriteDIASettingSerialization(String mzXMLFileName) {
        try {
            LogManager.getRootLogger().info("Writing DIA setting to file:" + FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_diasetting.ser...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(mzXMLFileName) + FilenameUtils.getBaseName(mzXMLFileName) + "_diasetting.ser", false);
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject(this);
            oos.close();
            fout.close();
        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }

    public static DIA_Setting ReadDIASettingSerialization(String filepath) {
        if (!new File(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser").exists()) {
            return null;
        }
        try {
            LogManager.getRootLogger().debug("Reading DIA setting from file:" + FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser...");

            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(filepath) + FilenameUtils.getBaseName(filepath) + "_diasetting.ser");
            ObjectInputStream in = new ObjectInputStream(fileIn);
            DIA_Setting setting = (DIA_Setting) in.readObject();
            in.close();
            fileIn.close();
            return setting;

        } catch (Exception ex) {
            LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;        
        }
    }

}
