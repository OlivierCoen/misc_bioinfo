import time
from urllib.request import urlretrieve
from goatools.obo_parser import GODag
from tempfile import gettempdir
from pathlib import Path
from unipressed import IdMappingClient, UniprotkbClient
from goatools.mapslim import mapslim
from collections import Counter
import seaborn as sns
import pandas as pd
from dataclasses import dataclass, field
from typing import ClassVar


@dataclass
class GOManager:

    GOTERM_CATEGORIES: ClassVar[list[str]] = ['molecular_function', 'biological_process', 'cellular_component']
    GO_OBO_PERMANENT_LINK: ClassVar[str] = 'http://purl.obolibrary.org/obo/go.obo'
    GOSLIM_GENERIC_OBO_PERMANENT_LINK: ClassVar[str] = 'https://current.geneontology.org/ontology/subsets/goslim_generic.obo'

    godag: GODag = field(init=False)
    goslim_dag: GODag = field(init=False)


    def __post_init__(self):
        obo_file = Path(gettempdir()) / 'go.obo'
        obo_slim_file = Path(gettempdir()) / 'goslim_generic.obo'
        print('Downloading obo file')
        urlretrieve(self.GO_OBO_PERMANENT_LINK, obo_file)
        print('Downloading goslim obo file')
        urlretrieve(self.GOSLIM_GENERIC_OBO_PERMANENT_LINK, obo_slim_file)
        self.godag = GODag(str(obo_file))
        self.goslim_dag = GODag(str(obo_slim_file))


    @staticmethod
    def get_uniprot_ids(tair_ids: list):
        request = IdMappingClient.submit(source="Gene_Name", dest="UniProtKB", ids=tair_ids)
        while request.get_status() != 'FINISHED':
            time.sleep(1)
        return list(request.each_result())


    @staticmethod
    def get_go_terms(uniprot_id: str):
        go_ids = []
        record_dict = UniprotkbClient.fetch_one(uniprot_id)
        if uniprotkb_crossrefs := record_dict['uniProtKBCrossReferences']:
            for cross_ref_dict in GOManager.get_specific_cross_references(uniprotkb_crossrefs, dbname='GO'):
                if 'id' in cross_ref_dict:
                    go_ids.append(cross_ref_dict['id'])
        return go_ids

    @staticmethod
    def get_specific_cross_references(uniprotkb_crossrefs: list, dbname: str):
        if not isinstance(uniprotkb_crossrefs, list):
            return []
        return [
            cross_ref_dict for cross_ref_dict in uniprotkb_crossrefs
            if cross_ref_dict.get('database') == dbname
        ]

    @staticmethod
    def merge_go_term_list(group: list):
        merged_go_term_list = []
        for group in group:
            merged_go_term_list += group
        return merged_go_term_list


    def get_goslims(self, go_term: str):
        return list(mapslim(go_term, self.godag, self.goslim_dag)[1])


    def get_goslim_mappings(self, row: pd.Series):
        goslim_lists = []
        for go_term in row['go_terms']:
            try:
                goslims = self.get_goslims(go_term)
            except:
                print(f'Could not get goslim from {go_term}')
            else:
                if goslims:
                    goslim_lists.append(goslims)
        goslim_name_lists = [[self.godag[goslim].name for goslim in goslim_list] for goslim_list in goslim_lists]
        return goslim_lists, goslim_name_lists


    def get_specific_goslim_mappings(self, row: pd.Series, namespace: str):
        filtered_goslims = list(set([
            goslim_list[-1] for goslim_list in row['goslim_terms']
            if self.godag[goslim_list[-1]].namespace == namespace
        ]))
        filtered_goslim_names = [self.godag[goslim].name for goslim in filtered_goslims]
        return filtered_goslims, filtered_goslim_names


    def go_term_corresponds_to_goslim(self, go_term: str, goslim_name: str):
        try:
            goslims = self.get_goslims(go_term)
        except:
            return False
        else:
            for goslim in goslims:
                if self.godag[goslim].name == goslim_name:
                    return True
            return False


    def keep_go_term_corresponding_to_goslim(self, go_terms: list, goslim_name: str):
        return [go_term for go_term in go_terms if self.go_term_corresponds_to_goslim(go_term, goslim_name)]


    def plot_counts(self, axes, namespace: str, goterms_lists: list, subplot_index: int):
        all_goterms = []
        for goterm_list in goterms_lists:
            all_goterms += goterm_list

        all_goterms = list(all_goterms)
        all_names = [self.godag[goslim].name for goslim in all_goterms]
        all_names_counter = Counter(all_names)

        namespace_df = pd.Series(dict(all_names_counter), name='counts').reset_index()
        namespace_df.columns = ['goslim', 'percentage']
        namespace_df['percentage'] = namespace_df['percentage'] / all_names_counter.total()
        namespace_df.sort_values(by='percentage', ascending=False, inplace=True)
        subdf = namespace_df.head(20)

        # plot
        palette = sns.color_palette("Accent", n_colors=len(subdf))
        bars = sns.barplot(x='goslim', y='percentage', hue='goslim', data=subdf, ax=axes[subplot_index], palette=palette, legend=False)
        for j, bar in enumerate(bars.patches):
            bar.set_color(palette[j])
        axes[subplot_index].set_title(namespace)
        axes[subplot_index].set_ylabel('Percentage')
        axes[subplot_index].set_xlabel(None)
        axes[subplot_index].set_xticklabels(axes[subplot_index].get_xticklabels(), rotation=45, ha='right')


    def plot_one_count(self, goterms_lists: list):
        all_goterms = []
        for goterm_list in goterms_lists:
            all_goterms += goterm_list

        all_goterms = list(all_goterms)
        all_names = [self.godag[goslim].name for goslim in all_goterms]
        all_names_counter = Counter(all_names)

        namespace_df = pd.Series(dict(all_names_counter), name='counts').reset_index()
        namespace_df.columns = ['goslim', 'percentage']
        namespace_df['percentage'] = namespace_df['percentage'] / all_names_counter.total()
        namespace_df.sort_values(by='percentage', ascending=False, inplace=True)
        subdf = namespace_df.head(20)

        # plot
        palette = sns.color_palette("Accent", n_colors=len(subdf))
        bars = sns.barplot(x='goslim', y='percentage', hue='goslim', data=subdf, palette=palette, legend=False)
        for j, bar in enumerate(bars.patches):
            bar.set_color(palette[j])