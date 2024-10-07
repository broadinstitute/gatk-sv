import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

function Feature({ Svg, title, description, colSize = 'col--4', contentAlignment = 'center' }) {
  const alignmentClass =
    contentAlignment === 'right' ? 'text--right' :
    contentAlignment === 'center' ? 'text--center' :
    'text--left';

  return (
    <div className={clsx('col', colSize)}>
      {Svg && (
        <div className="text--center">
          <Svg className={styles.featureSvg} role="img" />
        </div>
      )}
      <div className={clsx(alignmentClass, 'padding-horiz--md')}>
        {title && <h3>{title}</h3>}
        {description && <p>{description}</p>}
      </div>
    </div>
  );
}

export default Feature;
