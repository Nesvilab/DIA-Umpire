/* 
 * Author: Chih-Chiang Tsou <chihchiang.tsou@gmail.com>
 *             Nesvizhskii Lab, Department of Computational Medicine and Bioinformatics, 
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
            org.apache.log4j.Logger.getRootLogger().error(ExceptionUtils.getStackTrace(ex));
        }
    }
}
