// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

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
          // Please change this to your repo.
          // Remove this to remove the "edit this page" links.
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
              {
                label: 'Twitter',
                href: 'https://twitter.com/broadinstitute',
              },
            ],
          },
          {
            title: 'More',
            items: [
              {
                label: 'Talkowski lab',
                to: 'https://talkowski.mgh.harvard.edu',
              }
            ],
          },
        ],
        copyright: `Copyright © ${new Date().getFullYear()} BroadInstitute, Built with Docusaurus.`,
      },
      algolia: {
        appId: 'LI6UMHUDIS',
        apiKey: '97d929d265c25db1ed1816391a2a719a',
        indexName: 'gatk-sv',
        contextualSearch: true
    },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
