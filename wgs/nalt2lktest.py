


class LikelihoodTester(object):

	def __init__(self, region, groups):
		self.region = region
		self.groups = groups
		self.matrices = None


	def estimate(self, n_alt_data):

		lks = self.calculate_likelihood(n_alt_data)
		pass

	def evaluate(self):
		""" return the accuracy of the region """
		pass

	def prepare_matrix(self):
		""" create a set of matrices for each group """
		pass

	def calculate_likelihood(self, n_alt_data):
		""" calculate likelihood of each group """
		pass