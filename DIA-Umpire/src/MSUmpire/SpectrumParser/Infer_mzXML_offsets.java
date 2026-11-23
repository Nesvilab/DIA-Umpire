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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.eclipse.collections.impl.list.mutable.primitive.IntArrayList;
import org.eclipse.collections.impl.list.mutable.primitive.LongArrayList;

/**
 *
 * @author
 */
public class Infer_mzXML_offsets {

    /**
     * The scan elements of an mzXML file, in file order.
     */
    public static final class ScanOffsets {

        /**
         * Byte offset of each {@code <scan} tag, followed by the byte offset just past the last {@code </scan>}.
         */
        public final LongArrayList offsets = new LongArrayList();
        /**
         * The {@code num} attribute of each scan. It holds one entry fewer than {@link #offsets}.
         */
        public final IntArrayList scanNums = new IntArrayList();
    }

    private static final String SCAN_TAG = "<scan";
    // Attributes are separated by whitespace, so the leading \s keeps names such as "scanNum" from matching.
    private static final Pattern NUM_ATTRIBUTE = Pattern.compile("\\snum\\s*=\\s*([\"'])(\\d+)\\1");

    /**
     * Finds each scan tag and reads its scan number. Do not assume the n-th scan is scan n:
     * {@link mzXMLParser#to_mzXML} writes only spectra that have peaks, each under its scan number in the source file.
     * The file must use LF line endings, as the files written by to_mzXML do, because offsets are counted as line length + 1.
     */
    public static ScanOffsets infer_mzXML_scan_offsets(final InputStream is) throws IOException {
        final ScanOffsets ret = new ScanOffsets();
        long pos = 0;
        final BufferedReader bis = new BufferedReader(new InputStreamReader(is, StandardCharsets.ISO_8859_1));
        long posend = -1;
        while (true) {
            final String line = bis.readLine();
            if (line == null) {
                break;
            }
            int from = 0;
            int idx;
            while ((idx = IndexOfScanTag(line, from)) != -1) {
                // ParseElutionIndex reads the number right after "<scan num=", so the num attribute is on the line that opens the tag.
                final int tagEnd = line.indexOf('>', idx);
                ret.offsets.add(pos + idx);
                ret.scanNums.add(ParseScanNum(line.substring(idx, tagEnd == -1 ? line.length() : tagEnd), pos + idx));
                if (tagEnd == -1) {
                    break;
                }
                from = tagEnd + 1;
            }
            final int idxend = line.lastIndexOf("</scan>");
            if (idxend != -1) {
                posend = pos + idxend + "</scan>".length();
            }
            pos += line.length() + 1;
        }
        if (ret.offsets.isEmpty() || posend < ret.offsets.getLast()) {
            throw new RuntimeException("mzXML parse error");
        }
        ret.offsets.add(posend);
        return ret;
    }

    /**
     * Index of the next {@code <scan} start tag at or after {@code from}, skipping other elements whose names
     * start with "scan", such as {@code <scanOrigin}.
     */
    private static int IndexOfScanTag(final String line, int from) {
        while (true) {
            final int idx = line.indexOf(SCAN_TAG, from);
            if (idx == -1) {
                return -1;
            }
            final int next = idx + SCAN_TAG.length();
            if (next == line.length() || line.charAt(next) == '>' || Character.isWhitespace(line.charAt(next))) {
                return idx;
            }
            from = next;
        }
    }

    /**
     * Reads the num attribute from the text of a scan start tag.
     *
     * @param offset byte offset of the tag, for the error message
     */
    private static int ParseScanNum(final CharSequence scanTag, final long offset) throws IOException {
        final Matcher matcher = NUM_ATTRIBUTE.matcher(scanTag);
        if (!matcher.find()) {
            throw new IOException("No num attribute in the scan tag at byte offset " + offset + ": " + scanTag);
        }
        final long num;
        try {
            num = Long.parseLong(matcher.group(2));
        } catch (NumberFormatException ex) {
            throw new IOException("Scan number out of range at byte offset " + offset + ": " + matcher.group(2));
        }
        // Integer.MAX_VALUE is the key that marks the end of the last scan in the scan index.
        if (num >= Integer.MAX_VALUE) {
            throw new IOException("Scan number out of range at byte offset " + offset + ": " + num);
        }
        return (int) num;
    }
}
