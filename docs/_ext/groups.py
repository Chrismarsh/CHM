import re
import xml.etree.ElementTree as ET
import glob
from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst import Directive
from docutils.statemachine import ViewList
from sphinx.util.nodes import nested_parse_with_titles


class Provides(Directive):
    has_content = False
    required_arguments = 0
    final_argument_whitespace = False

    def should_remove(self, x):
        if 'filter' in x:
            return True
        if 'modules' in x:
            return True

        return False

    def run(self):
        env = self.state.document.settings.env  # sphinx.environment.BuildEnvironment
        config = env.config
        folder = config["groups_doxygen_dir"]

        files = glob.glob(folder + '/*.xml')
        files = [x for x in files if not self.should_remove(x)]


        # get list of modules
        try:
            filter_xml = ET.parse(folder+'/group__modules.xml' ).getroot()[0]
        except FileNotFoundError as e:
            return []

        # these are all of our filters/modules
        provides = {}
        for child in filter_xml:
            if child.tag == 'innerclass':
                refid = child.get('refid')

                mod = ET.parse(folder+'/%s.xml' % refid ).getroot()[0]

                file = '../' + mod.find('.//location').get('file')[:-4] + '.cpp' #assumes that the hpp and cpp file are side by side

                cppfile = open(file, 'r')
                filetext = cppfile.read()
                cppfile.close()
                matches = re.findall("""provides\(\"(.+)\"\)""", filetext)

                for m in matches:
                    if not m in provides:
                        provides[m] = []

                    provides[m].append(child.text) #these modules provide 'm'

        rst = ViewList()

        for k in sorted(provides):

            if len(provides[k]) == 0:
                continue

            rst.append(k, "", 0)
            rst.append('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', "", 0)

            for f in provides[k]:
                rst.append(f, "", 0)
                rst.append('+++++++++++++++++++++++++++++++++++++++++++', "", 0)
                rst.append('.. doxygenclass:: ' + f, "", 0)
                rst.append('\n' , "", 0)

        node = nodes.section()
        node.document = self.state.document
        nested_parse_with_titles(self.state, rst, node)

        return node.children

class Filter(Directive):
    has_content = False
    required_arguments = 0

    option_spec = {
        'group': directives.unchanged_required,
    }

    final_argument_whitespace = False


    def should_remove(self, x):
        if 'filter' in x:
            return True
        if 'modules' in x:
            return True

        return False

    def run(self):
        env = self.state.document.settings.env  # sphinx.environment.BuildEnvironment
        config = env.config
        folder = config["groups_doxygen_dir"]

        files = glob.glob(folder + '/group__*')
        files = [x for x in files if not self.should_remove(x)]

        what = self.options['group']

        # for our current group (filter or modules), find what is in it
        try:
            filter_xml = ET.parse(folder+'/group__%s.xml' % what).getroot()[0]
        except FileNotFoundError as e:
            return []

        # these are all of our filters/modules
        filters = set()
        for child in filter_xml:
            if child.tag == 'innerclass':
                filters.add(child.text)

        groups = {}
        groups['Uncategorized']=[]

        seen_filters =set()
        for f in files:
            xml = ET.parse(f).getroot()[0]

            title = ''

            tmp =[]
            for child in xml:
                if child.tag == 'innerclass' and child.text in filters:
                    tmp.append(child.text)
                    seen_filters.add(child.text)
                if child.tag == 'title':
                    title = child.text

            groups[title] = tmp

        groups['Uncategorized'] = [x for x in filters-seen_filters]

        rst = ViewList()

        for k in sorted(groups):

            if len(groups[k]) == 0:
                continue

            rst.append(k, "", 0)
            rst.append('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', "", 0)

            for f in groups[k]:
                rst.append(f, "", 0)
                rst.append('+++++++++++++++++++++++++++++++++++++++++++', "", 0)
                rst.append('.. doxygenclass:: ' + f, "", 0)
                rst.append('\n' , "", 0)

        node = nodes.section()
        node.document = self.state.document
        nested_parse_with_titles(self.state, rst, node)

        return node.children


def setup(app):

    app.add_directive("groups", Filter)
    app.add_directive("provides", Provides)
    app.add_config_value('groups_doxygen_dir', './doxygen/xml/', 'html')
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }