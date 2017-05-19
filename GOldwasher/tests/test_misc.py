from unittest import TestCase
from nose.tools import assert_raises
from nose.tools import raises

import pkg_resources
from .. import misc


class Test_read_annotation_file(TestCase):

    @classmethod
    def setUp(self):
        self.imp = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/id2desc.txt')
        self.obj = misc.Slicer(self.imp)
        self.out =  self.obj.read_annotation_file(self.imp)

    def test_is_dict(self):
        self.assertTrue(isinstance(self.out, dict))

    def test_is_not_zero(self):
        self.assertTrue(len(self.out) > 0)

    def test_IOError(self):
        assert_raises(IOError, self.obj.read_annotation_file, 'null.file')


class Test_writeoutAdditionalAnnotation(TestCase):

    @classmethod
    def setUp(self):
        self.imp = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/id2desc.txt')
        self.obj = misc.Slicer(self.imp)
        self.troll = pkg_resources.resource_filename('GOldwasher',
                                                   'forbidden/out.txt')
    @raises(IOError)
    def test_file_readable(self):
        self.obj.writeoutAdditionalAnnotation('null.file', '/tmp/nosetest.tmp')

    @raises(IOError)
    def test_file_writeable(self):
        self.obj.writeoutAdditionalAnnotation(self.imp, self.troll)


class Test_read_enrichment_tsv(TestCase):

    @classmethod
    def setUp(self):
        self.imp = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/goenr.txt')
        self.imp2 = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/keggenr.txt')
        self.obj = misc.Slicer(self.imp)

    def test_go_input(self):
        self.assertTrue(len(self.obj.read_enrichment_tsv(
                            self.imp ,'GO').columns) == 8)

    def test_kegg_input(self):
        self.assertTrue(len(self.obj.read_enrichment_tsv(
                            self.imp2 ,'KEGG Pathways').columns) == 7)


class Test_process_enrichment_values(TestCase):

    @classmethod
    def setUp(self):
        self.imp = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/id2desc.txt')
        self.pathdir = pkg_resources.resource_filename('GOldwasher',
                                                       'tests/input/')
        self.f = 'goenr.txt' 
        self.obj = misc.Slicer(self.imp)
        self.out = self.obj.process_enrichment_values(self.pathdir,
                                                           self.f,
                                                           0.01) 

    def test_is_dict(self):
        self.assertTrue(isinstance(self.out, dict))

    def test_output_length(self):
        self.assertTrue(len(self.out) > 0)


class Test_generate_support_js_file(TestCase):

    @classmethod
    def setUp(self):
        self.imp = pkg_resources.resource_filename('GOldwasher',
                                                   'tests/input/id2desc.txt')
        self.obj = misc.Slicer(self.imp)
        self.troll = pkg_resources.resource_filename('GOldwasher',
                                                   'forbidden/out.txt')
        self.f = 'lefile.txt'
        self.d = {"Phatr3_EG02093.t1": "Fake Functional Annotation 1", 
                  "Phatr3_J45821.t1": "Fake Functional Annotation 2", 
                  "Phatr3_J41413.t1": "Fake Functional Annotation 3",}

    @raises(IOError)
    def test_file_writeable(self):
        self.obj.generate_support_js_file(self.d, self.troll, self.f, 'mockvar')