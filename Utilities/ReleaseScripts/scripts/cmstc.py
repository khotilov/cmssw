"""CMS TagCollector Python API
"""

__author__ = "Miguel Ojeda"
__copyright__ = "Copyright 2010-2011, CERN CMS"
__credits__ = ["Miguel Ojeda"]
__license__ = "Unknown"
__maintainer__ = "Miguel Ojeda"
__email__ = "mojedasa@cern.ch"
__status__ = "Staging"

_tagcollector_url = 'https://cmstags.cern.ch/tc/'

import urllib
import urllib2
import cookielib
import json
import getpass

class TagCollector(object):
	"""CMS TagCollector Python API"""

	def __init__(self):
		self._url = _tagcollector_url
		self._cj = cookielib.CookieJar()
		self._opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(self._cj))

	def _open(self, page, params = None, data = None):
		url = self._url + page + '?'
		if params:
			url += urllib.urlencode(params)
		if data:
			data = urllib.urlencode(data)
		try:
			return self._opener.open(url, data).read()
		except urllib2.HTTPError as e:
			raise Exception(e.read().strip())

	def _openjson(self, page, params = None, data = None):
		return json.loads(self._open(page, params, data))

	def signIn(self, username, password):
		"""Sign in to TagCollector."""
		self._open('CmsTCLogin', data = {'username': username, 'password': password})

	def signInInteractive(self):
		"""Sign in to TagCollector, asking for the username and password."""
		username = raw_input('Username: ')
		password = getpass.getpass()
		self.signIn(username, password)

	def signOut(self):
		"""Sign out of TagCollector."""
		self._open('signOut')

	def getPackageTags(self, package):
		"""Get the tags published in TagCollector for a package.
		Note: TagCollector's published tags are a subset of CVS' tags."""
		return self._openjson('py_getPackageTags', {'package': package})

	def getPackageTagDescriptionFirstLine(self, package, tag):
		"""Get the first line of the descriptions of a tag."""
		return self._openjson('py_getPackageTagDescriptionFirstLine', {'package': package, 'tag': tag})

	def getPackageTagReleases(self, package, tag):
		"""Get the releases where a tag is."""
		return self._openjson('py_getPackageTagReleases', {'package': package, 'tag': tag})

	def getReleasesTags(self, releases, diff = False):
		"""Get the tags of one or more release.
		Optionally, return only the tags that differ between releases."""
		releases = json.dumps(releases)
		diff = json.dumps(diff)
		return self._openjson('py_getReleasesTags', {'releases': releases, 'diff': diff})

	def getReleaseTags(self, release):
		"""Get the tags of one release."""
		return self.getReleasesTags((release, ))

	def getTagsetTags(self, tagset):
		"""Get the tags of one tagset."""
		return self._openjson('py_getTagsetTags', {'tagset': tagset})

	def getTagsetInformation(self, tagset):
		"""Get the information of one tagset."""
		return self._openjson('py_getTagsetInformation', {'tagset': tagset})

	def getPendingApprovalTags(self, args):
		"""Prints Pending Approval tags of one or more releases,
                one or more tagsets, or both (i.e. it joins all the tags).
		Prints an error if several tags appear for a single package.
		Suitable for piping to addpkg (note: at the moment,
		addpkg does not read from stdin, use "-f" instead)."""
		args = json.dumps(args)
		return self._openjson('py_getPendingApprovalTags', {'args': args})

	def commentTagsets(self, tagset_ids, comment):
		"""Comment one or more tagsets.
		Requirement: Signed in."""
		tagset_ids = json.dumps(tagset_ids)
		if len(comment) < 1:
			raise Exception("Error: Expected a comment.")
		self._open('commentTagsets', {'tagset_ids': tagset_ids, 'comment': comment})

	def signTagsets(self, tagset_ids, comment = ''):
		"""Sign one or more tagsets.
		Requirement: Signed in as a L2."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('signTagsets', {'tagset_ids': tagset_ids, 'comment': comment})

	def signTagsetsAll(self, tagset_ids, comment = ''):
		"""Sign all one or more tagsets.
		Requirement: Signed in as a Release Manager."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('signTagsetsAll', {'tagset_ids': tagset_ids, 'comment': comment})

	def rejectTagsetsPendingSignatures(self, tagset_ids, comment = ''):
		"""Reject one or more tagsets Pending Signatures.
		Requirement: Signed in as a L2."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('rejectTagsetsPendingSignatures', {'tagset_ids': tagset_ids, 'comment': comment})

	def approveTagsets(self, tagset_ids, comment = ''):
		"""Approve one or more tagsets.
		Requirement: Signed in as a Release Manager for each tagset's release."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('approveTagsets', {'tagset_ids': tagset_ids, 'comment': comment})

	def bypassTagsets(self, tagset_ids, comment = ''):
		"""Bypass one or more tagsets.
		Requirement: Signed in as a Release Manager for each tagset's release."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('bypassTagsets', {'tagset_ids': tagset_ids, 'comment': comment})

	def rejectTagsetsPendingApproval(self, tagset_ids, comment = ''):
		"""Reject one or more tagsets Pending Approval.
		Requirement: Signed in as a Release Manager."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('rejectTagsetsPendingApproval', {'tagset_ids': tagset_ids, 'comment': comment})

	def removeTagsets(self, tagset_ids, comment = ''):
		"""Remove one or more tagsets from the History (i.e. stack of the release).
		Requirement: Signed in as a Release Manager."""
		tagset_ids = json.dumps(tagset_ids)
		self._open('removeTagsets', {'tagset_ids': tagset_ids, 'comment': comment})

	def getPackagesPendingApproval(self):
		"""Get New Package Requests which are Pending Approval."""
		return self._openjson('py_getPackagesPendingApproval')

	def getPackageManagersRequested(self, package):
		"""Get the Package Managers (administrators and developers) requested in a New Package Request."""
		return self._openjson('py_getPackageManagersRequested', {'package': package})

	def search(self, term):
		"""Searches for releases, packages, tagsets, users and categories.
		Requirement: Signed in."""
		return self._openjson('search', {'term': term})

	def approvePackage(self, package):
		"""Approve a New Package Request.
		Requirement: Signed in as a Creator (i.e. people in the top-level .admin/developers file).
		Warning: This does *not* create the package in CVS."""
		self._open('approveNewPackageRequest', {'package_name': package})

	def getIBs(self, filt = '', limit = 10):
		"""Get the name and creation date of Integration Builds.
		By default, it only returns the latest 10 IBs.
		Optionally, filter by name."""
		return self._openjson('py_getIBs', {'filt': filt, 'limit': limit})

