class FastqRead(object):
   def __init__(self, desc, seq, qual):
      self.desc = desc
      self.seq = seq
      self.qual = qual

   def trim(self, idx):
      if idx is None:
         return self
      else:
         return self.__class__(self.desc, self.seq[:idx], self.qual[:idx])

   def format_fastq(self):
      return "@{0}\n{1}\n+\n{2}\n".format(self.desc, self.seq, self.qual)

   @classmethod
   def parse(cls, f):
      for desc, seq, _, qual in _grouper(f, 4):
         desc = desc.rstrip()[1:]
         seq = seq.rstrip()
         qual = qual.rstrip()
         yield cls(desc, seq, qual)

def _grouper(iterable, n):
   "Collect data into fixed-length chunks or blocks"
   # grouper('ABCDEFG', 3) --> ABC DEF
   args = [iter(iterable)] * n
   return zip(*args)
