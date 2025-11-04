// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

import {themes as prismThemes} from 'prism-react-renderer';

/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'GATK-SV',
  tagline: 'A cloud-native pipeline for calling structural variations on short-read sequencing data',
  url: 'https://broadinstitute.github.io',
  baseUrl: '/gatk-sv/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'throw',
  //favicon: 'img/favicon.ico',

  // GitHub pages deployment config.
  organizationName: 'broadinstitute',
  projectName: 'gatk-sv',

  i18n: {
    defaultLocale: 'en',
    locales: ['en'],
  },

  presets: [
    [
      'classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      ({
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          editUrl:
            'https://github.com/broadinstitute/gatk-sv/tree/master/website',
        },
        /*blog: {
          showReadingTime: true,
          editUrl: '...',
        },*/
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
        gtag: {
          trackingID: 'G-9DQEYKHD1M',
          anonymizeIP: true,
        },
      }),
    ],
  ],

  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'GATK-SV',
        /*logo: {
          alt: 'GATK-SV logo',
          src: 'img/logo.svg',
        },*/
        items: [
          {
            type: 'doc',
            docId: 'intro',
            position: 'right',
            label: 'Documentation',
          },
          /*{
            to: '/blog',
            label: 'Blog',
            position: 'right'
          },*/
          {
			label: 'Questions',
			href: 'https://github.com/broadinstitute/gatk-sv/issues',
			position: 'right'
		  },
          {
            href: 'https://github.com/broadinstitute/gatk-sv',
            label: 'GitHub',
            position: 'right',
          },
        ],
      },
      footer: {
        style: 'dark',
        links: [
          {
            title: 'Docs',
            items: [
              {
                label: 'About',
                to: '/docs/intro',
              },
            ],
          },
          {
            title: 'Community',
            items: [
              {
                label: 'Github',
                href: 'https://github.com/broadinstitute/gatk-sv/discussions',
              },
            ],
          },
          {
            title: 'More',
            items: [
              {
                label: 'Talkowski lab',
                to: 'https://talkowski.mgh.harvard.edu',
              },
              {
                label: 'Broad Institute',
                to: 'https://www.broadinstitute.org'
              },
              {
                label: 'Center for Genomic Medicine',
                to: 'https://www.massgeneral.org/research/cgm'
              }
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} BroadInstitute, Built with Docusaurus.`,
      },
      algolia: {
        appId: 'FPHFCFLYVM',
        apiKey: '378661222bc7b372faaebf438f20eb5d',
        indexName: 'gatk-sv',
        contextualSearch: true
    },
    prism: {
        theme: prismThemes.github,
        darkTheme: prismThemes.vsDark,
        additionalLanguages: [
            'bash', 'powershell', 'json', 'json5', 'r', 'awk', 'jq', 'diff',
            'git', 'ignore', 'log', 'makefile', 'regex', 'docker']
      },
      docs: {
        sidebar: {
          hideable: true,
        }
      }
    }),

  themes: ['@docusaurus/theme-mermaid'],
  markdown: {
    mermaid: true,
  }
};

module.exports = config;
