class GTFExon:
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attributes):
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = self.parse_attributes(attributes) if isinstance(attributes, str) else attributes

    @staticmethod
    def parse_attributes(attr_string):
        """Parse the GTF attributes field into a dictionary."""
        attributes = {}
        fields = attr_string.strip().strip(';').split(';')
        for field in fields:
            field = field.strip()
            if not field:
                continue
            if '=' in field:
                key, value = field.split('=', 1)
            elif ' ' in field:
                key, value = field.split(' ', 1)
            else:
                key, value = field, ""
            attributes[key.strip()] = value.strip().strip('"').strip("'")
        return attributes
    
    def attr_str(self):
        """Return the attributes as a string."""
        return ' '.join([f'{k} "{v}";' for k, v in self.attributes.items()])

    def __str__(self):
        """Reconstruct the GTF line from the object's attributes."""
        attr_str = ' '.join([f'{k} "{v}";' for k, v in self.attributes.items()])
        return '\t'.join([
            self.seqname,
            self.source,
            self.feature,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.frame,
            attr_str
        ])

    @classmethod
    def from_gtf_line(cls, line):
        """Create a GTFExon object from a GTF line."""
        fields = line.strip().split('\t')
        if len(fields) != 9:
            raise ValueError("GTF line must have 9 fields")
        return cls(*fields)


class Transcript:
    def __init__(self, seqname, source, feature, start, end, score, strand, frame, attributes):
        self.seqname = seqname
        self.source = source
        self.feature = feature  # usually "transcript"
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = GTFExon.parse_attributes(attributes) if isinstance(attributes, str) else attributes
        self.exons = []

    def add_exon(self, exon):
        """Add a GTFExon object to the transcript."""
        if not isinstance(exon, GTFExon):
            raise TypeError("Exon must be an instance of GTFExon")
        self.exons.append(exon)

    def sort_exons(self):
        """Sort exons by start position, depending on strand."""
        self.exons.sort(key=lambda exon: exon.start)

    def exon_count(self):
        """Return number of exons."""
        return len(self.exons)
    
    def get_junctions(self):
        '''Get junctions from the exons.'''
        self.sort_exons()
        junctions = []
        if len(self.exons) >= 2:
            for i in range(len(self.exons) - 1):
                exon1 = self.exons[i]
                exon2 = self.exons[i + 1]
                junctions.append((exon1.end, exon2.start))
        return junctions
    
    def get_UJC(self):
        """Get UJC (Unique Junction chain) from the exons."""
        # UJC format: <chromosome>_<strand>_<junction1_start>_<junction1_end>_<junction2_start>_<junction2_end>_...
        chrom = self.seqname
        strand = self.strand
        junctions = self.get_junctions()
        if junctions:
            junctions_str = "_".join(map(lambda x: f"{x[0]}_{x[1]}", junctions))
        else:
            junctions_str = "mono-exon"
        return f"{chrom}_{strand}_{junctions_str}"

    def attr_str(self):
        """Return the attributes as a string."""
        return ' '.join([f'{k} "{v}";' for k, v in self.attributes.items()])
    
    def __str__(self):
        attr_str = ' '.join([f'{k} "{v}";' for k, v in self.attributes.items()])
        header = '\t'.join([
            self.seqname, self.source, self.feature,
            str(self.start), str(self.end),
            self.score, self.strand, self.frame,
            attr_str
        ])
        exon_lines = "\n".join(str(exon) for exon in self.exons)
        return f"{header}\n{exon_lines}\n"


def read_gtf_as_transcripts(gtf_path):
    transcripts = {}

    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # skip comments and empty lines

            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue  # malformed line

            feature = fields[2]
            attributes = GTFExon.parse_attributes(fields[8])

            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                continue  # skip lines without transcript_id

            if feature == "transcript":
                # Create and store the Transcript object
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = Transcript(*fields)
                else:
                    transcripts[transcript_id].start = min(transcripts[transcript_id].start, int(fields[3]))
                    transcripts[transcript_id].end = max(transcripts[transcript_id].end, int(fields[4]))
            elif feature == "exon":
                exon = GTFExon(*fields)
                if transcript_id not in transcripts:
                    # Auto-instantiate Transcript if no preceding "transcript" line exists
                    transcripts[transcript_id] = Transcript(
                        fields[0], fields[1], "transcript", fields[3], fields[4],
                        fields[5], fields[6], fields[7], fields[8]
                    )
                else:
                    transcripts[transcript_id].start = min(transcripts[transcript_id].start, int(fields[3]))
                    transcripts[transcript_id].end = max(transcripts[transcript_id].end, int(fields[4]))
                transcripts[transcript_id].add_exon(exon)

    for t in transcripts.values():
        t.sort_exons()
    return transcripts
