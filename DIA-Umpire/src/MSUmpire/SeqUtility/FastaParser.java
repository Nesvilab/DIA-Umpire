/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of computational medicine and bioinformatics, 
 *             University of Michigan, Ann Arbor
 *
 * Copyright 2014 University of Michigan, Ann Arbor, MI
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package MSUmpire.SeqUtility;

import MSUmpire.PSMDataStructure.EnzymeManager;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.nustaq.serialization.FSTObjectInput;
import org.nustaq.serialization.FSTObjectOutput;

/**
 *
 * @author Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 */
public class FastaParser implements Serializable{
    private static final long serialVersionUID = 19398249L;

    public HashMap<String, ProteinEntry> ProteinList;
    public HashMap<String, PeptideEntry> PeptideList;

    public void RemoveDecoy(String DecoyTag) {
        HashMap<String, ProteinEntry> newlist=new HashMap<>();
        for(ProteinEntry protein : ProteinList.values()){
            if(!(protein.ACC.startsWith(DecoyTag)|protein.ACC.endsWith(DecoyTag))){
                newlist.put(protein.ACC, protein);
            }
        }
        ProteinList=newlist;
    }
        
    public String GetProtSeq(String ProtID){        
        if(ProteinList.containsKey(ProtID)){
            return ProteinList.get(ProtID).Seq;
        }        
        for(String ID : ProteinList.keySet()){
            if(ID.contains(ProtID)){
                return ProteinList.get(ID).Seq;
            }
        }
        return null;
    }
    public class ProteinEntry implements Serializable{
        private static final long serialVersionUID = -2002064228475586294L;
        public String ACC;
        public String Des;
        public String Seq;
        public ArrayList<String> Peptides=new ArrayList<>();
    }
    
    public class PeptideEntry implements Serializable {
        private static final long serialVersionUID = -6343751134961266096L;
        public ArrayList<String> Proteins=new ArrayList<>();
        public String Sequence;
        public String Decoy;
    }
    
    public FastaParser(String filename){
        ProteinList=new HashMap<>();
        try {
            Parse(filename);
        } catch (IOException ex) {
           LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
    
    
     public static FastaParser FasterSerialzationRead(String Filename) throws FileNotFoundException {

        if (!new File(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename)+ ".FastaSer").exists()) {
            return null;
        }
        FastaParser fastareader=null;
        try {
            org.apache.logging.log4j.LogManager.getRootLogger().info("Loading fasta serialization to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename)+ ".FastaSer..");
            FileInputStream fileIn = new FileInputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename)+ ".FastaSer");
            FSTObjectInput in = new FSTObjectInput(fileIn);
            fastareader = (FastaParser) in.readObject();
            in.close();
            fileIn.close();
        } catch (Exception ex) {
            org.apache.logging.log4j.LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return null;
        }
        
        return fastareader;
    }

    public boolean FasterSerialzationWrite(String Filename) throws FileNotFoundException {
        try {
            org.apache.logging.log4j.LogManager.getRootLogger().info("Writing fasta serialization to file:" + FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename)+ ".FastaSer...");
            FileOutputStream fout = new FileOutputStream(FilenameUtils.getFullPath(Filename) + FilenameUtils.getBaseName(Filename)+ ".FastaSer", false);
            FSTObjectOutput out = new FSTObjectOutput(fout);
            out.writeObject(this);
            out.close();
            fout.close();
        } catch (Exception ex) {
            org.apache.logging.log4j.LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
            return false;
        }        
        return true;
    }
    
    private void Parse(String filename) throws FileNotFoundException, IOException {

        if(!new File(filename).exists()){
            org.apache.logging.log4j.LogManager.getRootLogger().warn("Fasta file cannot be found: "+filename);
        }
        BufferedReader reader = new BufferedReader(new FileReader(filename));
        String line = "";
        String ACC = "";
        String des = "";
        StringBuilder Seq = new StringBuilder();
        while ((line = reader.readLine()) != null) {
            if (line.startsWith(">")) {
                if (!"".equals(ACC) && !"".equals(Seq.toString())) {
                    ProteinEntry protein = new ProteinEntry();
                    protein.ACC = ACC;
                    protein.Des = des;
                    protein.Seq = Seq.toString().replace("*", "");
                    ProteinList.put(ACC, protein);
                    Seq = new StringBuilder();
                }
                ACC = line.trim().split(" ")[0].substring(1).split("\t")[0];
                des = line.replace(">" + ACC + " ", "");
            } else {
                if (!"".equals(line.trim())) {
                    Seq.append(line);
                }
            }
        }
        if (!"".equals(ACC) && !"".equals(Seq.toString())) {
            ProteinEntry protein = new ProteinEntry();
            protein.ACC = ACC;
            protein.Des = des;
            protein.Seq = Seq.toString().replace("*", "");
            ProteinList.put(ACC, protein);
        }
        reader.close();
    }
        
    public void digestion(int missedcleave, int minlength, int maxlength, String Decoytag) throws IOException {
        PeptideList=new HashMap<>();
        for (ProteinEntry protein : ProteinList.values()) {
            if(protein.ACC.startsWith(Decoytag)|protein.ACC.endsWith(Decoytag)){
                continue;
            }
            String Sequence = protein.Seq;            
            ArrayList<String> TheoPeptides = EnzymeManager.GetInstance().GetTrypsinNoP().digest(Sequence, missedcleave, minlength, maxlength);
            AddFirstMetDroppedPep(Sequence, missedcleave, minlength, maxlength, TheoPeptides);
            for (String pep : TheoPeptides) {
                if (!PeptideList.containsKey(pep)) {
                    PeptideEntry pepentry=new PeptideEntry();
                    pepentry.Sequence=pep;
                    String rev=new StringBuilder(pepentry.Sequence.subSequence(0, pepentry.Sequence.length()-1)).reverse().toString();
                    if (Sequence.indexOf(pepentry.Sequence) > 0) {
                        rev += String.valueOf(Sequence.charAt(Sequence.indexOf(pepentry.Sequence) - 1));
                    }
                    pepentry.Decoy=rev;
                    PeptideList.put(pep, pepentry);
                }
                PeptideList.get(pep).Proteins.add(protein.ACC);
                protein.Peptides.add(pep);
            }
        }
    }

    public void AddFirstMetDroppedPep(String Sequence, int missedcleave, int minlength, int maxlength, ArrayList<String> TheoPeptides) {
        if (String.valueOf(Sequence.charAt(0)).equals("M")) {
            int mc = 0;
            for (int i = 1; i < Sequence.length(); i++) {
                if (String.valueOf(Sequence.charAt(i)).equals("K") || String.valueOf(Sequence.charAt(i)).equals("R")) {
                    mc++;
                    if (mc > missedcleave) {
                        return;
                    }
                    String pep = Sequence.substring(1, i + 1);
                    if (pep.length() >= minlength && pep.length() <= maxlength && !TheoPeptides.contains(pep)) {
                        TheoPeptides.add(pep);
                    }
                }
            }
        }
    }
    
    
}
