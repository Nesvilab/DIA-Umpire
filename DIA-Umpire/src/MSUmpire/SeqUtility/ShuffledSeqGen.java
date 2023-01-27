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
package MSUmpire.SeqUtility;

import ExternalPackages.JAligner.Alignment;
import ExternalPackages.JAligner.NeedlemanWunschGotoh;
import ExternalPackages.JAligner.Sequence;
import ExternalPackages.JAligner.matrix.Matrix;
import ExternalPackages.JAligner.matrix.MatrixLoaderException;
import java.util.ArrayList;
import java.util.Collections;
import org.apache.commons.lang3.exception.ExceptionUtils;

/**
 * Generate shuffled sequence as decoy spectra for spectral library
 * Lam, H., Deutsch, E.W. & Aebersold, R. Artificial decoy spectral
libraries for false discovery rate estimation in spectral library
searching in proteomics. J. Proteome Res. 9, 605–610
(2010).
 * @author Chih-Chiang Tsou
 */
public class ShuffledSeqGen implements Runnable{
    
    public String seq;
    public String decoy;
    Matrix blosum62;
    public ShuffledSeqGen(String seq, Matrix blosum62){
        this.seq=seq;
        this.blosum62=blosum62;
    }
    public void Generate() throws MatrixLoaderException {
        //ProteinSequence s1 = new ProteinSequence(fragmentLib.Sequence);                
        Sequence s1 = new Sequence(seq);
        float similarity = 1f;
        int NoIterations = 10;        
        Sequence s2 = new Sequence(shuffle(seq));
        for (int i = 0; i < NoIterations; i++) {
            s2 = new Sequence(shuffle(s2.getSequence()));
            Alignment alignment = NeedlemanWunschGotoh.align(s1, s2, blosum62, 10f, 0.5f);
            similarity = (float) alignment.getSimilarity() / alignment.getSequence1().length;
            if (similarity < 0.7f) {
                decoy= s2.getSequence();
                return;
            } else if (i == NoIterations - 1) {                
                s2.setSequence(RandomSequenceGeneratorWoKPR.GetInstance().GetNext() + s2.getSequence());
                i = 0;
                if(s2.length()>s1.length()+3){
                    break;
                }
            }
        }
        decoy= shuffle(seq);
        return;
    }
    private String shuffle(String s) {
        ArrayList<Character> list = new ArrayList<>();
        ArrayList<Integer> KRP = new ArrayList<>();
        for (int i = 0; i < s.length(); i++) {
//            if(i>0 && s.charAt(i) == 'P' &&(s.charAt(i-1) == 'K' || s.charAt(i-1) == 'R')){
//                KRP.add(i);
//            }
            if (s.charAt(i) == 'K' || s.charAt(i) == 'R'|| s.charAt(i) == 'P') {
                KRP.add(i);
            } else {
                list.add(s.charAt(i));
            }
        }
        Collections.shuffle(list);
        String shuffledSeq = "";

        int offset = 0;
        for (int i = 0; i < s.length(); i++) {
            if (KRP.contains(i)) {
                shuffledSeq += String.valueOf(s.charAt(i));
                offset++;
            } else {
                shuffledSeq += String.valueOf(list.get(i - offset));
            }
        }
        return shuffledSeq;
    }

    @Override
    public void run() {
        try {
            Generate();
        } catch (MatrixLoaderException ex) {
            org.apache.logging.log4j.LogManager.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
