package org.ucb.c5.composition;

public class DNAPermutation implements Comparable<DNAPermutation>{
    private String seq;
    private double GC_content;
    private double hairpin;
    private int length;

    public DNAPermutation(String seq, double GC_content, double hairpin) {
        this.seq = seq;
        this.GC_content = GC_content;
        this.hairpin = hairpin;
    }

    public DNAPermutation(String seq, double GC_content) {
        this.seq = seq;
        this.GC_content = GC_content;
        this.length = seq.length();
    }

    public DNAPermutation(String seq) {
        this.seq = seq;
    }

    public String getSeq() {
        return seq;
    }

    public double getGC_content() {
        return GC_content;
    }

    public double getHairpin() {
        return hairpin;
    }

    public int getLength() {
        return length;
    }

    @Override
    public int compareTo(DNAPermutation p) {
        if (this.GC_content > p.getGC_content()) {
            return 1;
        }
        else if (p.getGC_content() > this.GC_content) {
            return 1;
        }
        return 0;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        DNAPermutation that = (DNAPermutation) o;

        if (Double.compare(that.GC_content, GC_content) != 0) return false;
        if (Double.compare(that.hairpin, hairpin) != 0) return false;
        if (length != that.length) return false;
        return seq != null ? seq.equals(that.seq) : that.seq == null;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = seq != null ? seq.hashCode() : 0;
        temp = Double.doubleToLongBits(GC_content);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(hairpin);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + length;
        return result;
    }

    @Override
    public String toString() {
        return "DNAPermutation{" +
                "seq='" + seq + '\'' +
                ", GC_content=" + GC_content +
                ", hairpin=" + hairpin +
                ", length=" + length +
                '}';
    }
}
