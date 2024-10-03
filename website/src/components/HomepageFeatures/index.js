import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';
import useBaseUrl from '@docusaurus/useBaseUrl';
import Link from '@docusaurus/Link';

/*
Note that imageUrl used in the following is optional.
*/

const cloudNativeFeatures = [
  {
    imageUrl: 'img/placeholder.svg',
  },
  {
    title: <>Cloud-Native</>,
    description: (
      <>
        Built from the ground up for large-scale analysis, optimized for
        cloud-native execution on the Google Cloud Platform.
      </>
    ),
    //imageUrl: 'img/....svg',
    buttonLink: "/docs/gs/runtime-env",
    buttonText: "Learn More"
  },
];

const scalableFeatures = [
  {
    title: (
      <>Built for Population-Scale Genomic Research</>
    ),
    description: (
      <>
        Used to produce high-quality SV call sets for large
        scale sequencing initiatives such as the Genome
        Aggregation Project (gnomAD)
      </>
    ),
    buttonLink: "/docs/gs/overview",
    buttonText: "Learn More"
  },
  {
    imageUrl: 'img/placeholder.svg',
  }
];

const accurateFeatures = [
  {
    imageUrl: 'img/placeholder.svg',
  },
  {
    title: <>Accurate</>,
    description: (
      <>
        Analyzes SV calls from multiple algorithms and evidence
        signatures to achieve high sensitivity and precision
      </>
    ),
  },
];

const accessibleFeatures = [
  {
    title: <>Accessible</>,
    description: (
      <>
        Leveraging the <a href="https://terra.bio" target="_blank" rel="noopener noreferrer">Terra platform </a>
        for scalable execution, secure environment, and seamless collaboration.
      </>
    ),
    buttons: [
      {
        buttonLink: 'https://app.terra.bio',
        buttonText: 'Joint Calling Workspace'
      },
      {
        buttonLink: 'https://app.terra.bio/#workspaces/help-gatk/GATK-Structural-Variants-Single-Sample',
        buttonText: 'Single-Sample Workspace'
      }
    ],
  },
  {
    imageUrl: 'img/terraLogo.png',
  },
];

const featuredProjects = {
  header: "Featured Projects",
  images: [
    {
      src: 'img/gnomad.png',
      link: 'https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants/',
      alt: 'gnomAD'
    },
    {
      src: 'img/allOfUs.png',
      link: 'https://support.researchallofus.org/hc/en-us/articles/27496716922900-All-of-Us-Short-Read-Structural-Variant-Quality-Report',
      alt: 'all-of-us'
    },
    {
      src: 'img/hgsvc.png',
      link: 'https://www.hgsvc.org/',
      alt: 'HGSVC'
    },
    {
      src: 'img/1kgp.png',
      link: 'https://www.internationalgenome.org/data-portal/data-collection/30x-grch38',
      alt: '1000-genomes-project'
    }
  ],
  //description: "..."
};

const featureList = [
  {
    title: 'Evidence Collection and QC',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        SV evidence is collected through multiple algorithms and passed through rigorous QC steps
        to filter technical outliers and improve call quality.
      </>
    ),
  },
  {
    title: 'SV Genotyping and Refinement',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        After SV discovery, the pipeline refines genotypes using a cohort-wide reference panel
        to improve the accuracy of detected SVs. It supports re-genotyping of
        copy-number variants (CNVs) for large-scale analyses.
      </>
    ),
  },
  {
    title: 'Downstream Annotation and Visualization',
    //Svg: require('@site/static/img/....svg').default,
    description: (
      <>
        The pipeline includes modules for downstream annotation and visualization,
        making it easier to interpret the results in the context of population and medical genetics.
      </>
    ),
  },
];

function WholeRowFeature({ imageUrl, title, description, buttons, contentAlignment, imageAlignment }) {
  const imgUrl = useBaseUrl(imageUrl);

  return (
    <div className={clsx('col col--6', styles.featureContainer, imageAlignment === 'right' ? styles.alignRight : styles.alignLeft)}>
      {/* Image Section */}
      {imgUrl && (
        <div className={styles.featureImage}>
          <img
            className={styles.largeFeatureImage}
            src={imgUrl}
            alt={title || ''}
          />
        </div>
      )}

      {/* Content Section */}
      <div className={clsx(styles.featureContent, contentAlignment === 'right' ? styles.alignRight : '')}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
        {buttons && buttons.length > 0 && (
          <div className={clsx(styles.buttonContainer, contentAlignment === 'right' ? styles.alignRight : '')}>
            {buttons.map((button, index) => (
              <Link key={index} className="button button--primary" to={button.buttonLink}>
                {button.buttonText}
              </Link>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}


function Feature({ Svg, title, description }) {
  return (
    <div className={clsx('col col--4')}>
      {Svg && (
        <div className="text--center">
          <Svg className={styles.featureSvg} role="img" />
        </div>
      )}
      <div className={clsx('text--center', 'padding-horiz--md')}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
      </div>
    </div>
  );
}

function FeatureGallery({ header, images, description }) {
  return (
    <div className={styles.featureGallery}>
      <h2 className={styles.header}>{header}</h2>

      <div className={styles.imageGallery}>
        {images.map((image, index) => (
          <a key={index} href={image.link} target="_blank" rel="noopener noreferrer">
            <img className={styles.galleryImage} src={image.src} alt={image.alt || `Image ${index + 1}`} />
          </a>
        ))}
      </div>

      {description && <p className={styles.description}>{description}</p>}
    </div>
  );
}

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
                  contentAlignment="right"
                  imageAlignment="left"
                />
              ))}
            </div>
          </div>
        </section>
      )}

      <section className={clsx(styles.features, styles.featuresAlt)}>
        <div className="container">
          <FeatureGallery {...featuredProjects} />
        </div>
      </section>

      {cloudNativeFeatures && cloudNativeFeatures.length > 0 && (
        <section className={styles.features}>
          <div className="container">
            <div className={clsx('row', 'single-feature-row')}>
              {cloudNativeFeatures.map((props, idx) => (
                <WholeRowFeature key={idx} {...props} />
              ))}
            </div>
          </div>
        </section>
      )}

      {scalableFeatures && scalableFeatures.length > 0 && (
        <section className={clsx(styles.features, styles.featuresAlt)}>
          <div className="container">
            <div className={clsx('row', 'single-feature-row')}>
              {scalableFeatures.map((props, idx) => (
                <WholeRowFeature key={idx} {...props} />
              ))}
            </div>
          </div>
        </section>
        )}

      {accurateFeatures && accurateFeatures.length > 0 && (
        <section className={clsx(styles.features)}>
          <div className="container">
            <div className={clsx('row', 'single-feature-row')}>
              {accurateFeatures.map((props, idx) => (
                <WholeRowFeature key={idx} {...props} />
              ))}
            </div>
          </div>
        </section>
        )}

      <section className={clsx(styles.features)}>
        <div className="container">
          <div className="row">
            {featureList.map((props, idx) => (
              <Feature key={idx} {...props} />
            ))}
          </div>
        </div>
      </section>
    </>
  );
}
