import operator
import sys
import parse_train
import math
import optparse


class HMM(object):
	def __init__(self, count_file, output_file):
		self._count_filename = count_file
		self._count_file = open(count_file,"r")
		self._test_file = open("gene.test","r")
		self._output_file = open(output_file,'w')
		self._emission = {}
		self._ngram = [{} for _ in range(3)]
		self._tag = None
		
	def parameters(self):
		for line in self._count_file:
			line_buffer = line.split()
			if line_buffer[1] == "WORDTAG":
				if line_buffer[3] not in self._emission:
					self._emission[line_buffer[3]] = {}
				self._emission[line_buffer[3]][line_buffer[2]] = float(line_buffer[0])

			if line_buffer[1].find("-GRAM") > 0:
				n = int(line_buffer[1][0])-1
				self._ngram[n][tuple(line_buffer[2:])] = float(line_buffer[0])


	def emission(self, x, y):
		if x not in self._emission:
			if self._count_filename.find("p3") == -1:
				x = "_RARE_"
			else:
				if parse_train.check_numeric(x):
					x = "_NUMERIC_"
				elif parse_train.check_all_capitals(x):
					x = "_ALL_CAPITALS_"
				elif parse_train.check_last_capital(x):
					x = "_LAST_CAPITAL_"
				else:
					x = "_RARE_"

		if y in self._emission[x]:
			return self._emission[x][y] / self._ngram[0][tuple([y])]
		else:
			return sys.float_info.min
	
	def transition(self, y1, y2, y3):
		if tuple([y1, y2, y3]) in self._ngram[2]:
			return self._ngram[2][tuple([y1, y2, y3])] / self._ngram[1][tuple([y1, y2])] 
		else:
			return sys.float_info.min

	def viterbi(self, sentence):
		n = len(sentence)
		K = ["I-GENE","O"]
		tag = ['']*(n+1)
		
		pi = {}
		back_pointer = {}

		
		pi[tuple([0, '*', '*'])] = 1

		for i in range(1,n+1):
			x = sentence[i-1]
			if i == 1:
				for v in K:
					pi[tuple([i, '*', v])] = pi[tuple([i-1, '*', '*'])]* self.transition('*', '*', v) * self.emission(x,v)
			elif i == 2:
				for u in K:
					for v in K:
						pi[tuple([i, u, v])] = pi[tuple([i-1, '*', u])]* self.transition('*', u, v) * self.emission(x,v)		
			else:
				for u in K:
					for v in K:
						likelihood = [pi[tuple([i-1, w, u])] * self.transition(w, u, v) * self.emission(x,v) for w in K]
						index, value = max(enumerate(likelihood), key=operator.itemgetter(1))
						pi[tuple([i, u, v])] = value
						back_pointer[tuple([i, u, v])] = K[index]

		likelihood = [pi[tuple([n, u, v])] * self.transition( u, v, 'STOP')  for u in K for v in K]
		index, value = max(enumerate(likelihood), key=operator.itemgetter(1))
		tag[n-1],tag[n] = K[index/2], K[index-(index/2)*2]

		for i in reversed(range(1,n-2+1)):
			tag[i] = back_pointer[tuple([i+2,tag[i+1],tag[i+2]])]
		del tag[0]
		self._tag = tag

	def forward_backward(self, sentence):
		n = len(sentence)
		K = ["I-GENE","O"]
		tag = ['']*(n+1)
		
		alpha = {}
		beta = {}

		#forward	
		alpha[tuple([0, '*', '*'])] = 1

		for i in range(1,n+1):
			x = sentence[i-1]
			if i == 1:
				for v in K:
					alpha[tuple([i, '*', v])] = alpha[tuple([i-1, '*', '*'])]* self.transition('*', '*', v) * self.emission(x,v)
			elif i == 2:
				for u in K:
					for v in K:
						alpha[tuple([i, u, v])] = alpha[tuple([i-1, '*', u])]* self.transition('*', u, v) * self.emission(x,v)		
			else:
				for u in K:
					for v in K:		
						alpha[tuple([i, u, v])] = sum([alpha[tuple([i-1, w, u])] * self.transition(w, u, v) * self.emission(x,v) for w in K])
						
		#backward
		for u in K:
			for v in K:
				beta[tuple([n, u, v])] = 1

		for i in reversed(range(1,n)):
			x = sentence[i]
			if i == 1:
				for v in K:
					beta[tuple([i, '*', v])] = sum( [beta[tuple([i+1, v, w])] * self.transition('*', v, w) * self.emission(x,w) for w in K] )
			else:			
				for u in K:
					for v in K:
						beta[tuple([i, u, v])] = sum( [beta[tuple([i+1, v, w])] * self.transition(u, v, w) * self.emission(x,w) for w in K] )

		#posterior
		#print alpha[tuple([i,'*', K[0]])], beta[tuple([i,'*', K[0]])]
		for i in range(1,n+1):
			if i == 1:
				if alpha[tuple([i,'*', K[0]])]*beta[tuple([i,'*', K[0]])] > alpha[tuple([i,'*', K[1]])] * beta[tuple([i,'*', K[1]])]:
					tag[i] = K[0]
				else:
					tag[i] = K[1]
			else:
				if sum([alpha[tuple([i,u, K[0]])]*beta[tuple([i,u, K[0]])] for u in K]) > sum([alpha[tuple([i,u, K[1]])]*beta[tuple([i,u, K[1]])] for u in K]):
					tag[i] = K[0]
				else:
					tag[i] = K[1]

		del tag[0]
		self._tag = tag

	def logadd(self, alpha):
		x = alpha[0]
		y = alpha[1]
		if y <= x:
			return x + math.log1p(math.exp(y-x))
		else:
			return y + math.log1p(math.exp(x-y))

	def forward_backward_logrithm(self, sentence):
		n = len(sentence)
		K = ["I-GENE","O"]
		tag = ['']*(n+1)
		
		alpha = {}
		beta = {}

		#forward	
		alpha[tuple([0, '*', '*'])] = 0

		for i in range(1,n+1):
			x = sentence[i-1]
			if i == 1:
				for v in K:
					alpha[tuple([i, '*', v])] = alpha[tuple([i-1, '*', '*'])] + math.log(self.transition('*', '*', v)) + math.log( self.emission(x,v))
			elif i == 2:
				for u in K:
					for v in K:
						alpha[tuple([i, u, v])] = alpha[tuple([i-1, '*', u])] + math.log( self.transition('*', u, v) ) + math.log( self.emission(x,v) )		
			else:
				for u in K:
					for v in K:
						summation = sum([ math.exp(alpha[tuple([i-1, w, u])]) * self.transition(w, u, v) * self.emission(x,v) for w in K])
						log_alpha = [ alpha[tuple([i-1, w, u])] + math.log(self.transition(w, u, v)) + math.log(self.emission(x,v)) for w in K ]
						alpha[tuple([i, u, v])] = self.logadd(log_alpha)
						
		#backward
		for u in K:
			for v in K:
				beta[tuple([n, u, v])] = 0

		for i in reversed(range(1,n)):
			x = sentence[i]
			if i == 1:
				for v in K:
					log_beta = [ beta[tuple([i+1, v, w])] + math.log(self.transition('*', v, w)) + math.log(self.emission(x,w)) for w in K ]
					beta[tuple([i, '*', v])] = self.logadd(log_beta)
			else:			
				for u in K:
					for v in K:
						log_beta = [ beta[tuple([i+1, v, w])] + math.log(self.transition(u, v, w)) + math.log(self.emission(x,w)) for w in K ]
						beta[tuple([i, u, v])] = self.logadd(log_beta)
		#posterior
		#print alpha[tuple([i,'*', K[0]])], beta[tuple([i,'*', K[0]])]
		for i in range(1,n+1):
			if i == 1:
				if alpha[tuple([i,'*', K[0]])]+beta[tuple([i,'*', K[0]])] > alpha[tuple([i,'*', K[1]])] + beta[tuple([i,'*', K[1]])]:
					tag[i] = K[0]
				else:
					tag[i] = K[1]
			else:
				former = [alpha[tuple([i,u, K[0]])]+beta[tuple([i,u, K[0]])] for u in K]
				later =  [alpha[tuple([i,u, K[1]])]+beta[tuple([i,u, K[1]])] for u in K]
				if self.logadd(former) >  self.logadd(later):
					tag[i] = K[0]
				else:
					tag[i] = K[1]

		del tag[0]
		self._tag = tag


	def gene_tagger_basline(self):
		output = self._output_file
		
		for line in self._test_file:
			if line == "\n":
				output.write("\n")
				continue
			line_buffer = line.split()
			word = line_buffer[0]

			if self.emission(word, "I-GENE") > self.emission(word, "O"):
				output.write(line_buffer[0]+" I-GENE\n")
			else:
				output.write(line_buffer[0]+" O\n")

	def gene_tagger_viterbi(self):
		output = self._output_file
		
		sentence = []
		for line in self._test_file:
			if line == "\n":
				self.viterbi(sentence)
				for i,t in enumerate(self._tag):
					output.write(sentence[i]+" "+t+"\n")
				output.write("\n")
				sentence = []
				continue
			line_buffer = line.split()
			sentence.append(line_buffer[0])

	def gene_tagger_forward_backward(self):
		output = self._output_file
		
		sentence = []
		for line in self._test_file:
			if line == "\n":
				self.forward_backward_logrithm(sentence)
				for i,t in enumerate(self._tag):
					output.write(sentence[i]+" "+t+"\n")
				output.write("\n")
				sentence = []
				continue
			line_buffer = line.split()
			sentence.append(line_buffer[0])

	
def get_options(args=None):
    optParser = optparse.OptionParser()
    optParser.add_option("-c", "--counts-file", dest="counts_file",
                            help="define the counts file (mandatory)")
    optParser.add_option("-o", "--output-file", dest="output_file",
                            help="define the output file (mandatory)")
    optParser.add_option("-a", type="int", dest="algorithm",
                         default=1, help="algorithm (default 1=baseline, 2=viterbi, 3=forward_backward)")
    
    (options, args) = optParser.parse_args(args=args)
    return options

if __name__ == "__main__":

	options = get_options()
	hmm = HMM(options.counts_file, options.output_file)
	hmm.parameters()
	if options.algorithm == 1:
		hmm.gene_tagger_basline()
	elif options.algorithm == 2:
		hmm.gene_tagger_viterbi()
	else:
		hmm.gene_tagger_forward_backward()
	


