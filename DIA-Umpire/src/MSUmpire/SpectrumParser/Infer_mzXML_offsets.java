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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;

/**
 *
 * @author
 */
public class Infer_mzXML_offsets {

    public static LongArrayList infer_mzXML_scanno(final InputStream is) throws IOException {
        final LongArrayList offsets = new LongArrayList();
        long pos = 0;
        final BufferedReader bis = new BufferedReader(new InputStreamReader(is, StandardCharsets.ISO_8859_1));
        long posend = -1;
        while (true) {
            final String line = bis.readLine();
            if (line == null) {
                break;
            }
            final int idx = line.indexOf("<scan");
            final int idxend = line.lastIndexOf("</scan>");
            if (idxend != -1) {
                posend = pos + idxend + "</scan>".length();
            }
            if (idx != -1) {
                offsets.add(pos + idx);
            }
            pos += line.length() + 1;
        }
        if (posend < offsets.getLast()) {
            throw new RuntimeException("mzXML parse error");
        }
        offsets.add(posend);
        return offsets;
    }
}
