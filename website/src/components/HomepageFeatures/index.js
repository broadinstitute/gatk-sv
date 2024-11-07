import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';
import useBaseUrl from '@docusaurus/useBaseUrl';
import Link from '@docusaurus/Link';
import FeatureGallery from './featureGallery.js';
import WholeRowFeature from './wholeRowFeature.js';
import Feature from './feature.js';

const accessibleFeatures = [
  {
    title: <>Terra integration</>,
    description: (
      <>
        Run your data with pre-configured workspaces on the secure
        <a href="https://terra.bio" target="_blank" rel="noopener noreferrer"> Terra platform </a>
      </>
    ),
    buttons: [
      {
        buttonLink: 'https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling',
        buttonText: 'Joint Calling Workspace'
      },
      {
        buttonLink: 'https://app.terra.bio/#workspaces/help-gatk/GATK-Structural-Variants-Single-Sample',
        buttonText: 'Single-Sample Workspace'
      }
    ],
    colSize: "col--12",
  }
];

const featuredProjects = {
  header: "Featured Projects",
  images: [
    {
      src: 'img/logo/gnomad.png',
      link: 'https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants/',
      alt: 'gnomAD'
    },
    {
      src: 'img/logo/allOfUs.png',
      link: 'https://support.researchallofus.org/hc/en-us/articles/27496716922900-All-of-Us-Short-Read-Structural-Variant-Quality-Report',
      alt: 'all-of-us'
    },
    {
      src: 'img/logo/hgsvc.png',
      link: 'https://www.hgsvc.org/',
      alt: 'HGSVC'
    },
    {
      src: 'img/logo/1kgp.png',
      link: 'https://www.internationalgenome.org/data-portal/data-collection/30x-grch38',
      alt: '1000-genomes-project'
    }
  ],
  //description: "..."
};

const organizations = {
  header: "Organizations",
  images: [
    {
      src: "img/logo/broad.png",
      link: "https://www.broadinstitute.org/",
      alt: "Broad Institute"
    },
    {
      src: "img/logo/talkowskiLab.png",
      link: "https://talkowski.mgh.harvard.edu/",
      alt: "Talkowski Lab"
    },
    {
      src: "img/logo/cgm.png",
      link: "https://cgm.massgeneral.org/",
      alt: "Center for Genomic Medicine"
    }
  ]
};

const characteristics = [
  {
    title: 'Population-scale capabilities',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        Used for SV discovery in flagship research studies including
        the Genome Aggregation Project (gnomAD) and All of Us.
      </>
    ),
  },
  {
    title: 'Sensitive and accurate',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        Ensemble calling with multiple SV discovery tools combined with
        joint genotyping maximize power, and ML-based variant adjudication filters poor quality variants.
      </>
    ),
  },
  {
    title: 'Cloud-native',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        Built for the Terra genomics platform, enabling scalability, collaboration,
        and reproducibility in a secure environment.
      </>
    ),
  },
];

export default function HomepageFeatures() {
  return (
    <>
      {accessibleFeatures && accessibleFeatures.length > 0 && (
        <section className={clsx(styles.features)}>
          <div className="container">
            <div className={clsx('row', 'single-feature-row')}>
              {accessibleFeatures.map((props, idx) => (
                <WholeRowFeature
                  key={idx}
                  {...props}
                  contentAlignment="center"
                  imageAlignment="center"
                />
              ))}
            </div>
          </div>
        </section>
      )}

      <section className={clsx(styles.featuresAlt)}>
        <div className="container">
          <div className="row">
            {characteristics.map((props, idx) => (
              <Feature key={idx} {...props} />
            ))}
          </div>
        </div>
      </section>

      <section className={clsx(styles.features, styles.features)}>
        <div className="container">
          <FeatureGallery {...featuredProjects} />
        </div>
      </section>

      <section className={clsx(styles.featuresAlt)}>
        <div className="container">
          <FeatureGallery {...organizations} />
        </div>
      </section>
    </>
  );
}
